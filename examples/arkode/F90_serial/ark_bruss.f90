!-----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
!-----------------------------------------------------------------
! LLNS/SMU Copyright Start
! Copyright (c) 2015, Southern Methodist University and 
! Lawrence Livermore National Security
!
! This work was performed under the auspices of the U.S. Department 
! of Energy by Southern Methodist University and Lawrence Livermore 
! National Laboratory under Contract DE-AC52-07NA27344.
! Produced at Southern Methodist University and the Lawrence 
! Livermore National Laboratory.
!
! All rights reserved.
! For details, see the LICENSE file.
! LLNS/SMU Copyright End
! Copyright (c) 2013, Southern Methodist University.
! All rights reserved.
! For details, see the LICENSE file.
!-----------------------------------------------------------------
! Example problem:
! 
! The following test simulates a brusselator problem from chemical 
! kinetics.  This is an ODE system with 3 components, Y = [u,v,w], 
! satisfying the equations,
!    du/dt = a - (w+1)*u + v*u^2
!    dv/dt = w*u - v*u^2
!    dw/dt = (b-w)/ep - w*u
! for t in the interval [0.0, 10.0], with initial conditions 
! Y0 = [u0,v0,w0].  We use the initial conditions and parameters 
!    u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
! Here, all three solution components exhibit a rapid transient 
! change during the first 0.2 time units, followed by a slow and 
! smooth evolution.
! 
! This program solves a the Fortran ODE test problem using the 
! FARKODE interface for the ARKode ODE solver module.
! 
! This program uses the IMEX ARK solver; here the 
! implicit systems are solved with a modified Newton iteration
! with the ARKDENSE dense linear solver.  The Jacobian routine 
! and right-hand side routines come from the file user-supplied 
! Jacobian routine.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
program driver
  ! Declarations
  implicit none

  ! general problem variables
  integer*8, parameter :: NEQ=3
  real*8,    parameter :: T0=0.d0, Tf=10.d0
  real*8    :: dTout, Tout, Tcur, rtol, atol, rout(6)
  integer   :: it, Nt, ier
  integer*8 :: iout(22)
  real*8, dimension(NEQ) :: y

  ! real/integer parameters to pass through to supplied functions
  !    ipar(1) -> unused
  !    rpar(1) -> "a" parameter
  !    rpar(2) -> "b" parameter 
  !    rpar(3) -> "ep" parameter
  integer*8 :: ipar
  real*8    :: rpar(3)

  ! solver parameters
  integer*8 :: order, adapt_method
  real*8    :: nlscoef

  !-----------------------
  ! set some solver parameters
  order = 3          ! 4th order method
  adapt_method = 0   ! PID-controller
  nlscoef = 1.d-2    ! Newton solver tolerance coefficient

  ! time-stepping information
  dTout = (Tf-T0)/10.d0    ! output time interval
  Nt = Tf/dTout + 0.5      ! number of outputs

  ! set initial conditions, problem parameters
  y(1) = 3.9d0     ! u0
  y(2) = 1.1d0     ! v0
  y(3) = 2.8d0     ! w0
  rpar(1) = 1.2d0  ! a
  rpar(2) = 2.5d0  ! b
  rpar(3) = 1.d-5  ! ep

  ! set tolerances according to problem specifications
  atol = 1.d-10
  rtol = 1.d-6
  
  ! initialize vector module
  call FNVInitS(4, NEQ, ier)
  if (ier < 0) then
     write(0,*) 'Error in FNVInitS = ',ier
     stop
  endif

  ! initialize ARKode solver to use IMEX integrator, scalar tolerances
  call FARKMalloc(T0, y, 2, 1, rtol, atol, &
                  iout, rout, ipar, rpar, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKMalloc = ',ier
     stop
  endif

  ! set integrator options
  call FARKSetIin('ORDER', order, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKSetIin = ',ier
     stop
  endif
  call FARKSetRin('NLCONV_COEF', nlscoef, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKSetIin = ',ier
     stop
  endif
  call FARKSetAdaptivityMethod(adapt_method, 1, 0, 0, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKSetAdaptMethod = ',ier
     stop
  endif

  ! specify use of dense linear solver
  call FARKDense(NEQ, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKDense = ',ier
     stop
  endif
  call FARKDenseSetJac(1, ier)
  if (ier < 0) then
     write(0,*) 'Error in FARKDenseSetJac = ',ier
     stop
  endif

  ! Open output stream for results, output comment line
  open(100, file='solution.txt')
  write(100,*) '# t u v w'

  ! output initial condition to disk 
  write(100,'(3x,4(es23.16,1x))') T0, y

  ! loop over time outputs
  Tout = T0
  Tcur = T0
  print *, '        t           u           v           w'
  print *, '  ----------------------------------------------------'
  print '(3x,4(es12.5,1x))', Tcur, y
  do it = 1,Nt

     Tout = min(Tout + dTout, Tf)           ! set next output time
     call FARKode(Tout, Tcur, y, 1, ier)    ! call solver
     if (ier < 0) then
        print *, 'Error at step ',it,', FARKode return flag =',ier
        exit
     end if

     ! output current solution
     print '(3x,4(es12.5,1x))', Tcur, y
     write(100,'(3x,4(es23.16,1x))') Tcur, y

  end do
  print *, '  ----------------------------------------------------'
  close(100)

  ! output solver statistics
  print *, '  '
  print *, 'Final Solver Statistics:'
  print '(2(A,i7),A)', '   Internal solver steps =', iout(3), &
       ' (attempted =', iout(6), ')'
  print '(2(A,i7))', '   Total RHS evals:  Fe =', iout(7), &
       ',  Fi =', iout(8)
  print '(A,i7)', '   Total linear solver setups =', iout(9)
  print '(A,i7)', '   Total RHS evals for setting up the linear system =', iout(17)
  print '(A,i7)', '   Total number of Jacobian evaluations =', iout(18)
  print '(A,i7)', '   Total number of Newton iterations =', iout(11)
  print '(A,i7)', '   Total number of nonlinear solver convergence failures =', iout(12)
  print '(A,i7)', '   Total number of error test failures =', iout(10)
  print *, '  '

  ! clean up
  call FARKFree()

end program driver
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Required subroutines for FARKODE interface
!-----------------------------------------------------------------

subroutine farkifun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! Implicit portion of the right-hand side of the ODE system
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,    intent(in)  :: t, rpar(3)
  integer*8, intent(in)  :: ipar(1)
  integer,   intent(out) :: ier
  real*8,    intent(in)  :: y(3)
  real*8,    intent(out) :: ydot(3)

  ! temporary variables
  real*8 :: u, v, w, a, b, ep

  ! set temporary values
  a  = rpar(1)
  b  = rpar(2)
  ep = rpar(3)
  u  = y(1)
  v  = y(2)
  w  = y(3)

  ! fill implicit RHS, set success flag
  ydot(1) = 0.d0
  ydot(2) = 0.d0
  ydot(3) = (b-w)/ep
  ier = 0
  
end subroutine farkifun
!-----------------------------------------------------------------

subroutine farkefun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! Explicit portion of the right-hand side of the ODE system
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,    intent(in)  :: t, rpar(3)
  integer*8, intent(in)  :: ipar(1)
  integer,   intent(out) :: ier
  real*8,    intent(in)  :: y(3)
  real*8,    intent(out) :: ydot(3)

  ! temporary variables
  real*8 :: u, v, w, a, b, ep

  ! set temporary values
  a  = rpar(1)
  b  = rpar(2)
  ep = rpar(3)
  u  = y(1)
  v  = y(2)
  w  = y(3)

  ! fill explicit RHS, set success flag
  ydot(1) = a - (w+1.d0)*u + v*u*u
  ydot(2) = w*u - v*u*u
  ydot(3) = -w*u
  ier = 0
  
end subroutine farkefun
!-----------------------------------------------------------------

subroutine farkdjac(neq,t,y,fy,DJac,h,ipar,rpar,wk1,wk2,wk3,ier)
!-----------------------------------------------------------------
! Jacobian computation routine
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,    intent(in)  :: t, h, rpar(3)
  integer*8, intent(in)  :: neq, ipar(1)
  integer,   intent(out) :: ier
  real*8,    intent(in), dimension(neq) :: y, fy, wk1, wk2, wk3
  real*8,    intent(out) :: DJac(neq,neq)

  ! temporary variables
  real*8 :: u, v, w, a, b, ep

  ! set temporary values
  a  = rpar(1)
  b  = rpar(2)
  ep = rpar(3)
  u  = y(1)
  v  = y(2)
  w  = y(3)

  ! fill implicit Jacobian, set success flag
  DJac = 0.d0
  DJac(3,3) = -1.d0/ep
  ier = 0
    
end subroutine farkdjac
!-----------------------------------------------------------------
