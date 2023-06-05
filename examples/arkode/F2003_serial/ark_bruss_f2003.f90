! ------------------------------------------------------------------
! Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
!                modified by Jean M. Sexton @ LBL
!                modified by Daniel M. Margolis @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2023, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
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
! with the SUNDENSE linear solver.  The Jacobian routine and
! right-hand side routines are supplied.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! number of equations
  integer(c_long), parameter :: neq = 3

  ! ODE parameters
  double precision, parameter, dimension(neq) :: y0 = (/ 3.9d0, 1.1d0, 2.8d0 /)
  double precision, parameter :: a = 1.2d0
  double precision, parameter :: b = 2.5d0
  double precision, parameter :: ep = 1.d-5

contains

  ! ----------------------------------------------------------------
  ! ExpRhsFn provides the right hand side explicit function for the
  ! ODE: dy1/dt = f1(t,y1,y2,y3)
  !      dy2/dt = f2(t,y1,y2,y3)
  !      dy3/dt = f3(t,y1,y2,y3)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function ExpRhsFn(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn                  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr) :: user_data                     ! user-defined data

    ! local data
    real(c_double) :: u, v, w

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yvec(:)
    real(c_double), pointer, dimension(neq) :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! set temporary values
    u  = yvec(1)
    v  = yvec(2)
    w  = yvec(3)

    ! fill RHS vector
    fvec(1) = a - (w + 1.d0)*u + v*u*u
    fvec(2) = w*u - v*u*u
    fvec(3) = -w*u

    ! return success
    ierr = 0
    return

  end function ExpRhsFn


  ! ----------------------------------------------------------------
  ! ImpRhsFn provides the right hand side implicit function for the
  ! ODE: dy1/dt = f1(t,y1,y2,y3)
  !      dy2/dt = f2(t,y1,y2,y3)
  !      dy3/dt = f3(t,y1,y2,y3)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function ImpRhsFn(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn                  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr) :: user_data                     ! user-defined data

    ! local data
    real(c_double) :: u, v, w

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: yvec(:)
    real(c_double), pointer, dimension(neq) :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! set temporary values
    u  = yvec(1)
    v  = yvec(2)
    w  = yvec(3)

    ! fill RHS vector
    fvec(1) = 0.d0
    fvec(2) = 0.d0
    fvec(3) = (b-w)/ep

    ! return success
    ierr = 0
    return

  end function ImpRhsFn

  ! ----------------------------------------------------------------
  ! Jac: The Jacobian function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function Jac(t, sunvec_y, sunvec_f, sunmat_J, user_data, &
       sunvec_t1, sunvec_t2, sunvec_t3) result(ierr) bind(C,name='Jac')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr),    value :: user_data ! user-defined data
    type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! pointers to data in SUNDIALS vector and matrix
    real(c_double), pointer :: J(:,:)


    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    J(1:3, 1:3) => FSUNDenseMatrix_Data(sunmat_J)

    ! fill Jacobian entries
    J(1:3, 1:3) = 0.d0
    J(3,3) = -1.d0/ep

    ! return success
    ierr = 0
    return

  end function Jac

end module ode_mod

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKode module
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep module
  use fsundials_nvector_mod      ! Fortran interface to the generic N_Vector
  use fsundials_matrix_mod       ! Fortran interface to the generic SUNMatrix
  use fsundials_linearsolver_mod ! Fortran interface to the generic SUNLinearSolver
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsundials_context_mod      ! Fortran interface to SUNContext
  use ode_mod                    ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)    :: ctx                      ! SUNDIALS context for the simulation
  real(c_double) :: tstart                   ! initial time
  real(c_double) :: tend                     ! final time
  real(c_double) :: rtol, atol               ! relative and absolute tolerance
  real(c_double) :: dtout                    ! output time interval
  real(c_double) :: tout                     ! output time
  real(c_double) :: tcur(1)                  ! current time
  integer(c_int) :: imethod, idefault, pq    ! time step adaptivity parameters
  real(c_double) :: adapt_params(3)          ! time step adaptivity parameters
  integer(c_int) :: ierr                     ! error flag from C functions
  integer(c_int) :: nout                     ! number of outputs
  integer(c_int) :: outstep                  ! output loop counter
  integer(c_long):: mxsteps                  ! max num steps

  real(c_double), parameter :: nlscoef = 1.d-2  ! non-linear solver coefficient
  integer(c_int), parameter :: order = 3        ! method order

  type(N_Vector),        pointer :: sunvec_y    ! sundials vector
  type(SUNMatrix),       pointer :: sunmat_A    ! sundials matrix
  type(SUNLinearSolver), pointer :: sunls       ! sundials linear solver
  type(c_ptr)                    :: arkode_mem  ! ARKODE memory
  real(c_double),        pointer :: yvec(:)     ! underlying vector

  !======= Internals ============

  ! create the SUNDIALS context
  ierr = FSUNContext_Create(c_null_ptr, ctx)

  ! initialize ODE
  tstart = 0.0d0
  tend   = 10.0d0
  tcur   = tstart
  tout   = tstart
  dtout  = 1.0d0
  nout   = ceiling(tend/dtout)

  ! create SUNDIALS N_Vector
  sunvec_y => FN_VNew_Serial(neq, ctx)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  yvec => FN_VGetArrayPointer(sunvec_y)

  ! initialize solution vector
  yvec = y0

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_funloc(ExpRhsFn), c_funloc(ImpRhsFn), tstart, sunvec_y, ctx)
  if (.not. c_associated(arkode_mem)) print *,'ERROR: arkode_mem = NULL'

  ! Tell ARKODE to use a dense linear solver.
  sunmat_A => FSUNDenseMatrix(neq, neq, ctx)
  if (.not. associated(sunmat_A)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  sunls => FSUNLinSol_Dense(sunvec_y, sunmat_A, ctx)
  if (.not. associated(sunls)) then
     print *, 'ERROR: sunls = NULL'
     stop 1
  end if

  ierr = FARKStepSetLinearSolver(arkode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    print *, 'Error in FARKStepSetLinearSolver'
    stop 1
  end if

  ierr = FARKStepSetJacFn(arkode_mem, c_funloc(Jac))
  if (ierr /= 0) then
    print *, 'Error in FARKStepSetJacFn'
    stop 1
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-10

  ierr = FARKStepSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
     print *, 'Error in FARKStepSStolerances, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepSetOrder(arkode_mem, order)
  if (ierr /= 0) then
    print *, 'Error in FARKStepSetOrder'
    stop 1
  end if

  ierr = FARKStepSetNonlinConvCoef(arkode_mem, nlscoef)
  if (ierr /= 0) then
    print *, 'Error in FARKStepSetNonlinConvCoef'
    stop 1
  end if

!!$  mxsteps = 5000
!!$  ierr = FARKStepSetMaxNumSteps(arkode_mem, mxsteps)
!!$  if (ierr /= 0) then
!!$    print *, 'Error in FARKStepSetNonlinConvCoef'
!!$    stop 1
!!$  end if

  imethod = 0
  idefault = 1
  pq = 0
  adapt_params = 0.d0
  ierr = FARKStepSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params)
  if (ierr /= 0) then
     print *, 'Error in FARKStepSetAdaptivityMethod, ierr = ', ierr, '; halting'
     stop 1
  end if

    ! Open output stream for results, output comment line
  open(100, file='solution.txt')
  write(100,*) '# t u v w'

  ! output initial condition to disk
  write(100,'(3x,4(es23.16,1x))') tstart, yvec


  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '       t           u           v           w       '
  print *, ' ---------------------------------------------------'
  print '(2x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)
  do outstep = 1,nout

     ! call ARKStep
     tout = min(tout + dtout, tend)
     ierr = FARKStepEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
     if (ierr /= 0) then
        print *, 'Error in FARKStepEvolve, ierr = ', ierr, '; halting'
        stop 1
     endif

     ! output current solution
     print '(2x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)
     write(100,'(3x,4(es23.16,1x))') tcur, yvec(1), yvec(2), yvec(3)

  enddo
  print *, ' ----------------------------------------------------'
  close(100)

  ! diagnostics output
  call ARKStepStats(arkode_mem)

  ! clean up
  call FARKStepFree(arkode_mem)
  call FN_VDestroy(sunvec_y)
  call FSUNMatDestroy(sunmat_A)
  ierr = FSUNLinSolFree(sunls)
  ierr = FSUNContext_Free(ctx)

end program main


! ----------------------------------------------------------------
! ARKStepStats
!
! Print ARKODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine ARKStepStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: ierr          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nst_a(1)      ! num steps attempted
  integer(c_long) :: nfe(1)        ! num explicit function evals
  integer(c_long) :: nfi(1)        ! num implicit function evals
  integer(c_long) :: nlinsetups(1) ! num linear solver setups
  integer(c_long) :: netfails(1)   ! num error test fails

  real(c_double)  :: hinused(1)    ! initial step size
  real(c_double)  :: hlast(1)      ! last step size
  real(c_double)  :: hcur(1)       ! step size for next step
  real(c_double)  :: tcur(1)       ! internal time reached

  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations

  !======= Internals ============

  ierr = FARKStepGetNumSteps(arkode_mem, nsteps)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumSteps, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumStepAttempts, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumRhsEvals, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetActualInitStep(arkode_mem, hinused)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetActualInitStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetLastStep(arkode_mem, hlast)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetLastStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetCurrentStep(arkode_mem, hcur)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetCurrentStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetCurrentTime(arkode_mem, tcur)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetCurrentTime, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumLinSolvSetups(arkode_mem, nlinsetups)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumLinSolvSetups, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumErrTestFails(arkode_mem, netfails)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumErrTestFails, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumNonlinSolvIters(arkode_mem, nniters)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumNonlinSolvIters, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumNonlinSolvConvFails, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumJacEvals(arkode_mem, njacevals)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumJacEvals, retval = ', ierr, '; halting'
     stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)'    ,'Total internal steps taken    =',nsteps
  print '(4x,A,i9)'    ,'Total internal steps attempts =',nst_a
  print '(4x,A,i9)'    ,'Total rhs exp function call   =',nfe
  print '(4x,A,i9)'    ,'Total rhs imp function call   =',nfi
  print '(4x,A,i9)'    ,'Num lin solver setup calls    =',nlinsetups
  print '(4x,A,i9)'    ,'Num error test failures       =',netfails
  print '(4x,A,es12.5)','First internal step size      =',hinused
  print '(4x,A,es12.5)','Last internal step size       =',hlast
  print '(4x,A,es12.5)','Next internal step size       =',hcur
  print '(4x,A,es12.5)','Current internal time         =',tcur
  print '(4x,A,i9)'    ,'Num nonlinear solver iters    =',nniters
  print '(4x,A,i9)'    ,'Num nonlinear solver fails    =',nncfails
  print *, ' '

  return

end subroutine ARKStepStats
