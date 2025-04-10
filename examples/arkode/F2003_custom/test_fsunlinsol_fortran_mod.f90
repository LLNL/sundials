! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2025, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! Program to test custom fsunlinsol_fortran_mod implementation
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Utility module for error-checking
! ------------------------------------------------------------------
module fsunlinsol_test_mod
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fsunlinsol_fortran_mod
  use fsunmatrix_fortran_mod
  use fnvector_fortran_mod
  implicit none

contains
  ! ------------------------------------------------------------------
  integer(c_int) function check_vector(sunvec_x, sunvec_y, tol, Nvar, N) result(failure)

    implicit none
    real(c_double), value :: tol
    integer(c_int64_t), value :: Nvar, N
    Type(N_Vector) :: sunvec_x, sunvec_y
    Type(FVec), pointer :: x, y
    integer(c_int64_t) :: i, j

    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)
    failure = 0
    do j = 1, N
      do i = 1, Nvar
        if (dabs(x%data(i, j) - y%data(i, j)) > tol) then
          failure = 1
        end if
      end do
    end do

    if (failure == 1) then
      print *, '  '
      print *, 'check_vector failure, differences:'
      print *, '   blk    idx       x        y       diff'
      print *, '  --------------------------------------------'
      do j = 1, N
        do i = 1, Nvar
          if (dabs(x%data(i, j) - y%data(i, j)) > tol) then
            print '(2x,2(i4,3x),3(es9.2,1x))', j, i, x%data(i, j), &
              y%data(i, j), dabs(x%data(i, j) - y%data(i, j))
          end if
        end do
      end do
      print *, '  --------------------------------------------'
      print *, '  '
    end if

  end function check_vector

end module fsunlinsol_test_mod

! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fsunlinsol_test_mod

  !======= Declarations =========
  implicit none

  ! local variables
  type(c_ptr)                      :: sunctx
  integer(c_int)                   :: fails, retval, j, k
  integer(c_int64_t), parameter :: N = 1000
  integer(c_int64_t), parameter :: Nvar = 50
  type(SUNMatrix), pointer   :: sA
  type(FMat), pointer   :: A
  type(SUNLinearSolver), pointer   :: LS
  type(FLinSol), pointer   :: S
  type(N_Vector), pointer   :: sX, sY, sB
  type(FVec), pointer   :: X, B

  !======= Internals ============

  ! initialize failure total
  fails = 0

  ! create SUNDIALS context
  fails = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! create new matrices and vectors
  sX => FN_VNew_Fortran(Nvar, N, sunctx)
  if (.not. associated(sX)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  X => FN_VGetFVec(sX)

  sY => FN_VNew_Fortran(Nvar, N, sunctx)
  if (.not. associated(sY)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sB => FN_VNew_Fortran(Nvar, N, sunctx)
  if (.not. associated(sB)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  B => FN_VGetFVec(sB)

  sA => FSUNMatNew_Fortran(Nvar, N, sunctx)
  if (.not. associated(sA)) then
    print *, 'ERROR: sunmat = NULL'
    stop 1
  end if
  A => FSUNMatGetFMat(sA)

  ! fill A and X with uniformly-distributed random numbers in [0,1)
  call random_number(X%data)
  call random_number(A%data)

  ! update A to scale by 1/Nvar, and 1 to anti-diagonal of each diagonal block
  do k = 1, N
    A%data(:, :, k) = A%data(:, :, k)/Nvar
    do j = 1, Nvar
      A%data(Nvar - j + 1, j, k) = A%data(Nvar - j + 1, j, k) + 1.d0
    end do
  end do

  ! compute B = A*X
  retval = FSUNMatMatvec(sA, sX, sB)
  if (retval /= SUN_SUCCESS) then
    print *, 'ERROR: FSUNMatMatvec fail'
    stop 1
  end if

  ! create custom linear solver
  LS => FSUNLinSolNew_Fortran(Nvar, N, sunctx)
  if (.not. associated(LS)) then
    print *, 'ERROR: sunlinsol = NULL'
    stop 1
  end if
  S => FSUNLinSolGetFLinSol(LS)

  ! test SUNLinSolGetType
  if (FSUNLinSolGetType(LS) /= SUNLINEARSOLVER_DIRECT) then
    fails = fails + 1
    print *, '>>> FAILED test -- FSUNLinSolGetType'
    print *, '    Unrecognized vector type', FSUNLinSolGetType(LS)
  else
    print *, 'PASSED test -- FSUNLinSolGetType'
  end if

  ! test SUNLinSolSetup
  retval = FSUNLinSolSetup(LS, sA)
  if (retval /= SUN_SUCCESS) then
    fails = fails + 1
    print *, '>>> FAILED test -- FSUNLinSolSetup'
  else
    print *, 'PASSED test -- FSUNLinSolSetup'
  end if

  ! test SUNLinSolSolve
  call FN_VConst(0.d0, sY)
  retval = FSUNLinSolSolve(LS, sA, sY, sB, 1.d-9)
  if ((check_vector(sX, sY, 1.d-15*Nvar*Nvar, Nvar, N) /= 0) &
      .or. (retval /= SUN_SUCCESS)) then
    fails = fails + 1
    print *, '>>> FAILED test -- FSUNLinSolSolve'
  else
    print *, 'PASSED test -- FSUNLinSolSolve'
  end if

  ! free solver, matrix and vectors
  call FSUNMatDestroy(sA)
  call FN_VDestroy(sX)
  call FN_VDestroy(sY)
  call FN_VDestroy(sB)
  retval = FSUNLinSolFree(LS)
  if (retval /= 0) then
    fails = fails + 1
    print *, '>>> FAILED test -- FSUNLinSolFree'
  else
    print *, 'PASSED test -- FSUNLinSolFree'
  end if

  ! free SUNDIALS context
  retval = FSUNContext_Free(sunctx)

  ! print results
  if (fails > 0) then
    print '(a,i3,a)', 'FAIL: FSUNLinSol module failed ', fails, ' tests'
    stop 1
  else
    print *, 'SUCCESS: FSUNLinSol module passed all tests'
  end if
  print *, '  '

end program main
