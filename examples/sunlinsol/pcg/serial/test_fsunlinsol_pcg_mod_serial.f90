! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! This file tests the Fortran 2003 interface to the SUNDIALS
! PCG SUNLinearSolver implementation. We run three tests to
! exercise the solver:
!   1. simple tridiagonal system (no preconditioning)
!   2. simple tridiagonal system (jacobi preconditioning)
!   3. tridiagonal system w/ scale vector s (no preconditioning)
! -----------------------------------------------------------------

module test_fsunlinsol_pcg_serial
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  implicit none

  integer(C_LONG), private, parameter :: N = 100
  integer(C_INT),  private, parameter :: pretype = 1     ! Preconditioning type (1 or 2)
  integer(C_INT),  private, parameter :: maxl    = 500   ! maxium Krylov subspace dimension (> 0)
  real(C_DOUBLE),  private, parameter :: tol     = 1e-13 ! solver tolerance

  type, private :: UserData
    integer(C_LONG) :: N
    type(N_Vector), pointer  :: d, s
  end type

contains

  integer(C_INT) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod
    use fnvector_serial_mod
    use fsunlinsol_pcg_mod
    use test_sunlinsol

    implicit none

    type(SUNLinearSolver), pointer :: LS         ! test linear solver
    type(SUNMatrix),       pointer :: A          ! dummy SUNMatrix (set to null)
    type(N_Vector),        pointer :: x, xhat, b ! test vectors
    type(N_Vector),        pointer :: s2         ! dummy scaling vector (set to null)
    type(UserData),        pointer :: probdata   ! problem data
    real(C_DOUBLE),        pointer :: xdata(:)   ! x vector data
    real(C_DOUBLE)                 :: tmpr       ! temporary real value
    integer(C_LONG)                :: j
    integer(C_INT)                 :: tmp

    ! setup
    fails = 0

    A  => null()
    s2 => null()

    x    => FN_VNew_Serial(N)
    xhat => FN_VNew_Serial(N)
    b    => FN_VNew_Serial(N)

    allocate(probdata)
    probdata%N = N
    probdata%d => FN_VNew_Serial(N)
    probdata%s => FN_VNew_Serial(N)

    ! fill xhat vector with uniform random data in [1, 2)
    xdata => FN_VGetArrayPointer(xhat)
    do j=1, N
      call random_number(tmpr)
      xdata(j) = ONE + tmpr
    end do

    ! fill Jacobi vector with matrix diagonal
    call FN_VConst(FIVE, probdata%d)

    ! create PCG linear solver
    LS => FSUNLinSol_PCG(x, pretype, maxl)

    ! run initialization tests
    fails = fails + Test_FSUNLinSolGetType(LS, SUNLINEARSOLVER_ITERATIVE, 0)
    fails = fails + Test_FSUNLinSolSetATimes(LS, c_loc(probdata),&
                                             c_funloc(ATimes), 0)
    fails = fails + Test_FSUNLinSolSetPreconditioner(LS,&
                                                     c_loc(probdata),&
                                                     c_funloc(PSetup),&
                                                     c_funloc(PSolve),&
                                                     0)
    fails = fails + Test_FSUNLinSolSetScalingVectors(LS, probdata%s,&
                                                     s2, 0)
    fails = fails + Test_FSUNLinSolInitialize(LS, 0)
    fails = fails + Test_FSUNLinSolSpace(LS, 0)

    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_PCG module, initialization'
    else
      print *, 'SUCCESS: FSUNLinSol_PCG module, initialization tests passed'
      print *, ''
    end if

    ! Test 1: simple Poisson-like solve (no preconditioning)

    ! set scaling vectors
    call FN_VConst(ONE, probdata%s)

    ! fill x vector with scaled version
    call FN_VDiv(xhat, probdata%s, x)

    ! fill b vector with result of matrix-vector product
    fails = ATimes(c_loc(probdata), x, b)
    if (fails /= 0) then
      print *, 'FAIL ATimes, problem 1'
      return
    end if

    ! Run tests with this setup
    fails = fails + FSUNLinSol_PCGSetPrecType(LS, PREC_NONE);
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0);
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0);
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0);
    fails = fails + Test_FSUNLinSolNumIters(LS, 0);
    fails = fails + Test_FSUNLinSolResNorm(LS, 0);
    fails = fails + Test_FSUNLinSolResid(LS, 0);

    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_PCG module, problem 1'
    else
      print *, 'SUCCESS: FSUNLinSol_PCG module, problem 1, all tests passed'
      print *, ''
    end if

    ! Test 2: simple Poisson-like solve (Jacobi preconditioning)

    ! set scaling vectors
    call FN_VConst(ONE, probdata%s)

    ! fill x vector with scaled version
    call FN_VDiv(xhat, probdata%s, x)

    ! fill b vector with result of matrix-vector product
    fails = ATimes(c_loc(probdata), x, b)
    if (fails /= 0) then
      print *, 'FAIL ATimes, problem 2'
      return
    end if

    ! Run tests with this setup
    fails = fails + FSUNLinSol_PCGSetPrecType(LS, pretype);
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0);
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0);
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0);
    fails = fails + Test_FSUNLinSolNumIters(LS, 0);
    fails = fails + Test_FSUNLinSolResNorm(LS, 0);
    fails = fails + Test_FSUNLinSolResid(LS, 0);

    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_PCG module, problem 2'
    else
      print *, 'SUCCESS: FSUNLinSol_PCG module, problem 2, all tests passed'
      print *, ''
    end if

    ! Test 3: Poisson-like solve w/ scaled rows (no preconditioning)

    ! set scaling vectors
    xdata => FN_VGetArrayPointer(probdata%s)
    do j=1, N
      call random_number(tmpr)
      xdata(j) = ONE + 1000.0d0*tmpr
    end do

    ! fill x vector with scaled version
    call FN_VProd(xhat, probdata%s, x)

    ! fill b vector with result of matrix-vector product
    fails = ATimes(c_loc(probdata), x, b)
    if (fails /= 0) then
      print *, 'FAIL ATimes, problem 3'
      return
    end if

    ! Run tests with this setup
    fails = fails + FSUNLinSol_PCGSetPrecType(LS, PREC_NONE);
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0);
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, tol, 0);
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0);
    fails = fails + Test_FSUNLinSolNumIters(LS, 0);
    fails = fails + Test_FSUNLinSolResNorm(LS, 0);
    fails = fails + Test_FSUNLinSolResid(LS, 0);

    if (fails /= 0) then
      print *, 'FAIL: FSUNLinSol_PCG module, problem 3'
    else
      print *, 'SUCCESS: FSUNLinSol_PCG module, problem 3, all tests passed'
      print *, ''
    end if

    ! cleanup
    tmp = FSUNLinSolFree(LS)
    call FN_VDestroy(x)
    call FN_VDestroy(xhat)
    call FN_VDestroy(b)
    call FN_VDestroy(probdata%d)
    call FN_VDestroy(probdata%s)

  end function unit_tests

  integer(C_INT) function ATimes(udata, vvec, zvec) result(ret) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use test_utilities

    implicit none

    type(C_PTR), value :: udata
    type(N_Vector) :: vvec, zvec
    type(UserData), pointer :: probdata
    real(C_DOUBLE), pointer :: v(:), z(:), s(:)
    integer(C_LONG) :: i, N

    call c_f_pointer(udata, probdata)

    v => FN_VGetArrayPointer(vvec)
    z => FN_VGetArrayPointer(zvec)
    s => FN_VGetArrayPointer(probdata%s)
    N = probdata%N

    ! perform product at the left domain boundary (note: v is zero at the boundary)
    z(1) = (FIVE*v(1)/s(1) - v(2)/s(2))/s(1)

    ! iterate through interior of local domain, performing product
    do i = 2, N-1
      z(i) = (-v(i-1)/s(i-1) + FIVE*v(i)/s(i) - v(i+1)/s(i+1))/s(i)
    end do

    ! perform product at the right domain boundary (note: v is zero at the boundary)
    z(N) = (-v(N-1)/s(N-1) + FIVE*v(N)/s(N))/s(N)

    ret = 0
  end function ATimes

  integer(C_INT) function PSetup(udata) result(ret) bind(C)
    use, intrinsic :: iso_c_binding
    type(C_PTR), value :: udata
    ret = 0
  end function PSetup

  integer(C_INT) function PSolve(udata, rvec, zvec, tol, lr) &
      result(ret) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use test_utilities

    implicit none

    type(C_PTR),    value   :: udata
    type(N_Vector)          :: rvec, zvec
    real(C_DOUBLE)          :: tol
    integer(C_INT)          :: lr
    type(UserData), pointer :: probdata
    real(C_DOUBLE), pointer :: r(:), z(:), d(:), s(:)
    integer(C_LONG)         :: i, N

    call c_f_pointer(udata, probdata)

    r => FN_VGetArrayPointer(rvec)
    z => FN_VGetArrayPointer(zvec)
    d => FN_VGetArrayPointer(probdata%d)
    s => FN_VGetArrayPointer(probdata%s)
    N = probdata%N

    do i=1, N
      z(i) = s(i) * s(i) * r(i) / d(i)
    end do

    ret = 0
  end function PSolve

end module

integer(C_INT) function check_vector(X, Y, tol) result(failure)
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_pcg_serial
  use test_utilities

  implicit none
  type(N_Vector)  :: x, y
  real(C_DOUBLE)  :: tol, maxerr
  integer(C_LONG) :: i, xlen, ylen
  real(C_DOUBLE), pointer :: xdata(:), ydata(:)

  failure = 0

  xdata => FN_VGetArrayPointer(x)
  ydata => FN_VGetArrayPointer(y)

  xlen = FN_VGetLength(x)
  ylen = FN_VGetLength(y)

  if (xlen /= ylen) then
    print *, 'FAIL: check_vector: different data array lengths'
    failure = 1
    return
  end if

  do i = 1, xlen
    failure = failure + FNEQTOL(xdata(i), ydata(i), FIVE*tol*abs(xdata(i)))
  end do

  if (failure > 0) then
    maxerr = ZERO
    do i = 1, xlen
      maxerr = max(abs(xdata(i)-ydata(i))/abs(xdata(i)), maxerr)
    end do
    write(*,'(A,E14.7,A,E14.7,A)') &
      "FAIL: check_vector failure: maxerr = ", maxerr, "  (tol = ", FIVE*tol, ")"
  end if

end function check_vector

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_pcg_serial

  !======== Declarations ========
  implicit none
  integer(C_INT) :: fails = 0

  !============== Introduction =============
  print *, 'PCG SUNLinearSolver Fortran 2003 interface test'
  print *, ''

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: ', fails, ' unit tests failed'
    stop 1
  else
    print *,'SUCCESS: all unit tests passed'
  end if
end program main
