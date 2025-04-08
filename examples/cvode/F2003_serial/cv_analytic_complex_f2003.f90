! ------------------------------------------------------------------
! Programmer(s): Mustafa Aggul @ SMU and Cody Balos @ LLNL
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
! This program solves a complex-valued ODE system using CVODE with
! the BDF method and a Newton iteration with the SPTFQMR linear solver.
!
! The ODE system is:
!   du/dt = (t - 1)*(t - w)*i
!   dv/dt = 2*u*i - v + 4
!   dw/dt = v + t^2*exp(-t)*i - exp(-t)*i + 1
!
! The true solution is:
!   u = -t*exp(-t) + 2*i
!   v = -t^2*exp(-t)*i
!   w = exp(-t)*i + t
!
! ------------------------------------------------------------------

module cv_analytic_complex_mod

  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod

  implicit none

  integer(c_int64_t), parameter :: neq = 3

contains

  ! ----------------------------------------------------------------
  ! f: The CVODE RHS operator function
  ! ----------------------------------------------------------------
  integer(c_int) function f(t, sunvec_y, sunvec_f, user_data) &
    result(retval) bind(C, name='f')

    use, intrinsic :: iso_c_binding

    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: sunvec_y
    type(N_Vector)        :: sunvec_f
    type(c_ptr), value    :: user_data

    complex(c_double_complex), pointer, dimension(neq) :: yval(:)
    complex(c_double_complex), pointer, dimension(neq) :: fval(:)

    yval => FN_VGetArrayPointer(sunvec_y)
    fval => FN_VGetArrayPointer(sunvec_f)

    fval(1) = (t - 1.0d0)*(t - yval(3))*(0.0d0, 1.0d0)
    fval(2) = 2.0d0*yval(1)*(0.0d0, 1.0d0) - yval(2) + 4.0d0
    fval(3) = yval(2) + t**2*exp(-t)*(0.0d0, 1.0d0) - exp(-t)*(0.0d0, 1.0d0) + 1.0d0

    retval = 0
    return

  end function f
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Solution: Compute the true solution
  ! ----------------------------------------------------------------
  integer(c_int) function Solution(t, sunvec_u) &
    result(retval) bind(C, name='Solution')

    use, intrinsic :: iso_c_binding

    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: sunvec_u

    complex(c_double_complex), pointer, dimension(neq) :: uval(:)

    uval => FN_VGetArrayPointer(sunvec_u)

    uval(1) = -t*exp(-t) + (0.0d0, 2.0d0)
    uval(2) = -t**2*exp(-t)*(0.0d0, 1.0d0)
    uval(3) = exp(-t)*(0.0d0, 1.0d0) + t

    retval = 0
    return

  end function Solution
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! SolutionError: Compute the solution error
  ! ----------------------------------------------------------------
  integer(c_int) function SolutionError(t, sunvec_u, sunvec_e) &
    result(retval) bind(C, name='SolutionError')

    use, intrinsic :: iso_c_binding

    implicit none

    real(c_double), value :: t
    type(N_Vector)        :: sunvec_u
    type(N_Vector)        :: sunvec_e

    real(c_double) :: error_norm
    complex(c_double_complex), pointer, dimension(neq) :: uval(:)
    complex(c_double_complex), pointer, dimension(neq) :: eval(:)

    uval => FN_VGetArrayPointer(sunvec_u)
    eval => FN_VGetArrayPointer(sunvec_e)

    retval = Solution(t, sunvec_e)
    call FN_VLinearSum((1.0d0, 0.0d0), sunvec_u, (-1.0d0, 0.0d0), sunvec_e, sunvec_e)
    error_norm = FN_VMaxNorm(sunvec_e)
    write(*, '(A, ES12.5)') "    Max-norm of the error is ", error_norm

    retval = 0
    return

  end function SolutionError
  ! ----------------------------------------------------------------

end module cv_analytic_complex_mod
! ------------------------------------------------------------------

program main

  use, intrinsic :: iso_c_binding
  use fcvode_mod
  use fnvector_serial_mod
  use fsunlinsol_spgmr_mod
  use cv_analytic_complex_mod

  implicit none

  real(c_double) :: dTout, t0, tf, tout, tret(1)
  real(c_double) :: reltol, abstol
  integer(c_int) :: retval, iout, nt

  type(N_Vector), pointer :: sunvec_y
  type(N_Vector), pointer :: sunvec_true
  type(N_Vector), pointer :: sunvec_err
  type(SUNLinearSolver), pointer :: sunlinsol_LS
  type(SUNMatrix), pointer:: sunmat_A
  type(c_ptr) :: cvode_mem, sunctx, stdout_file

  complex(c_double_complex), pointer, dimension(neq) :: yval(:)

  ! Problem parameters
  t0 = 0.0d0
  tf = 5.0d0
  tout = 1.0d0
  dTout = 1.0d0
  nt = ceiling(tf / dTout)
  reltol = 1.0d-6
  abstol = 1.0d-10

  ! Create SUNDIALS context
  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! Create C file pointer for stdout
  retval = FSUNDIALSFileOpen("stdout", "w+", stdout_file)

  ! Allocate vectors
  sunvec_y => FN_VNew_Serial(neq, sunctx)
  sunvec_true => FN_VClone(sunvec_y)
  sunvec_err => FN_VClone(sunvec_y)

  yval => FN_VGetArrayPointer(sunvec_y)
  yval(1) = -t0*exp(-t0) + (0.0d0, 2.0d0)
  yval(2) = -t0**2*exp(-t0)*(0.0d0, 1.0d0)
  yval(3) = exp(-t0)*(0.0d0, 1.0d0) + t0

  print *, ""
  print *, "Analytic ODE test problem:"
  print '(A, ES12.5, A, ES12.5)', "    reltol = ", reltol, ",  abstol = ", abstol
  print *, ""

  ! Create CVODE memory
  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  retval = FCVodeInit(cvode_mem, c_funloc(f), t0, sunvec_y)
  retval = FCVodeSStolerances(cvode_mem, reltol, abstol)

  ! Create linear solver
  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_y, SUN_PREC_NONE, 10, sunctx)
  sunmat_A => null()
  retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)

  ! Main time-stepping loop
  print *, "     t               u                      v                      w"
  print *, "   ----------------------------------------------------------------------------"
  print '(f8.3, " | ", 3(1x, f8.5, " + ", f8.5, "i", 1x))', tret(1), &
    real(yval(1)), aimag(yval(1)), real(yval(2)), aimag(yval(2)), &
    real(yval(3)), aimag(yval(3))
  do iout = 0, nt-1
    retval = FCVode(cvode_mem, tout, sunvec_y, tret, CV_NORMAL)
    print '(f8.3, " | ", 3(1x, f8.5, " + ", f8.5, "i", 1x))', tret(1), &
      real(yval(1)), aimag(yval(1)), real(yval(2)), aimag(yval(2)), &
      real(yval(3)), aimag(yval(3))
    tout = tout + 1.0d0
  end do
  print *, "   ----------------------------------------------------------------------------"
  retval = SolutionError(tf, sunvec_y, sunvec_err)
  print *, "   ----------------------------------------------------------------------------"
  print *, ""

  retval = FCVodePrintAllStats(cvode_mem, stdout_file, SUN_OUTPUTFORMAT_TABLE)

  ! Free memory
  call FCVodeFree(cvode_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_true)
  call FN_VDestroy(sunvec_err)
  retval = FSUNDIALSFileClose(stdout_file)
  retval = FSUNContext_Free(sunctx)

end program main
