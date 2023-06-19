! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
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
! Simple diagonal test with Fortran interface, using user-supplied
! preconditioner setup and solve routines (supplied in Fortran).
!
! This example does a basic test of the solver by solving the
! system:
!    f(u) = 0  for
!    f(u) = u(i)^2 - i^2
!
! No scaling is done.
! An approximate diagonal preconditioner is used.
!
! ------------------------------------------------------------------

module diag_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  integer(c_int),  parameter :: probsize = 128 

  integer(c_int)  :: ierr, i
  integer(c_long) :: iout(16)
  real(c_double)  :: p(probsize), rout(2), u(probsize), pscale(probsize), constr(probsize)
  integer(c_int),  parameter :: globalstrat = 0
  integer(c_int),  parameter :: prectype = 2 
  integer(c_int),  parameter :: maxl = 10
  integer(c_int),  parameter :: maxlrst = 2
  integer(c_long), parameter :: neq = probsize
  integer(c_long), parameter :: msbpre = 5
  real(c_double),  parameter :: fnormtol = 1.0d-5
  real(c_double),  parameter :: scsteptol = 1.0d-4

contains

  ! ----------------------------------------------------------------
  ! func: The nonlinear residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function func(sunvec_u, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_u  ! solution N_Vector
    type(N_Vector)       :: sunvec_f  ! LHS N_Vector
    type(c_ptr),   value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: uu(:), ff(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
    ff(1:neq) => FN_VGetArrayPointer(sunvec_f)

    ! loop over domain, computing our system f(u) = 0
    do i = 1,neq

       ! applying the constraint f(u) = u(i)^2 - i^2
       ff(i) = uu(i)*uu(i) - i*i
    end do


    ! return success
    ierr = 0
    return

  end function func

  integer(c_int) function kpsetup(sunvec_u, sunvec_s, sunvec_f, &
       sunvec_fs, user_data) result(ierr) bind(C)
    
    !======= Inclusions ===========
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_u  ! solution N_Vector
    type(N_Vector)       :: sunvec_s  ! scaling N_Vector
    type(N_Vector)       :: sunvec_f  ! LHS N_Vector
    type(N_Vector)       :: sunvec_fs ! LHS scaling N_Vector
    type(c_ptr),   value :: user_data ! user-defined data

    integer(c_int) :: i, ierr

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: udata(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    udata(1:neq) => FN_VGetArrayPointer(sunvec_u)

    ! loop over domain
    do i = 1,neq

      ! setup preconditioner
      p(i) = 0.5d0 / (udata(i) + 5.0d0)
    end do


    ! return success
    ierr = 0
    return

  end function kpsetup

  integer(c_int) function kpsolve(sunvec_u, sunvec_s, sunvec_f, &
       sunvec_fs, sunvec_v, user_data) result(ierr) bind(C)
    
    !======= Inclusions ===========
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_u  ! solution N_Vector
    type(N_Vector)       :: sunvec_s  ! scaling N_Vector
    type(N_Vector)       :: sunvec_f  ! LHS N_Vector
    type(N_Vector)       :: sunvec_fs ! LHS scaling N_Vector
    type(N_Vector)       :: sunvec_v  ! LHS scaling N_Vector
    type(c_ptr),   value :: user_data ! user-defined data

    integer(c_int) :: i, ierr

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(neq) :: udata(:)
    real(c_double), pointer, dimension(neq) :: v(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    udata(1:neq) => FN_VGetArrayPointer(sunvec_u)
    v(1:neq)     => FN_VGetArrayPointer(sunvec_v)

    ! loop over domain
    do i = 1,neq

      ! preconditioner solver
      v(i) = v(i) * p(i)
    end do


    ! return success
    ierr = 0
    return

  end function kpsolve

end module diag_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_context_mod
  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to SPGMR SUNLinearSolver
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_matrix_mod       ! Fortran interface to generic SUNmatrix
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use diag_mod                   ! problem-defining functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: ftol, fnorm(1)
  integer(c_int)  :: ierr

  type(c_ptr)                    :: sunctx        ! sundials context
  type(N_Vector),        pointer :: sunvec_u      ! sundials vectors
  type(N_Vector),        pointer :: sunvec_s      ! sundials vectors
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver

  type(c_ptr) :: kmem ! KINSOL memory

  real(c_double), dimension(nx,ny) :: u, scale

  !======= Internals ============

  ! -------------------------
  ! Print problem description

  print *, " "
  print *, "Example program fkinDiagon_kry:"
  print *, "   This FKINSOL example solves a 128 eqn diagonal algebraic system."
  print *, " Its purpose is to demonstrate the use of the Fortran interface in"
  print *, " in a serial environment."
  print *, " "
  print *, "Solution method: KIN_none"
  print '(a,i3)', "Problem size: neq = ", neq

  ! -------------------------
  ierr = FSUNContext_Create(c_null_ptr, sunctx)
  if (ierr /= 0) then
    print *, 'ERROR in FSUNContext_Create'
    stop 1
  end if

  ! -------------------------
  ! Set initial guess, and disable scaling

  do i = 1,neq
    u(i) = 2.0d0
    pscale(i) = 1.0d0
    constr(i) = 0.0d0
  continue

  ! -------------------------
  ! Create vectors for solution and scales

  sunvec_u => FN_VMake_Serial(neq, u, sunctx)
  if (.not. associated(sunvec_u)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_s => FN_VMake_Serial(neq, pscale, sunctx)
  if (.not. associated(sunvec_u)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  ! -------------------------
  ! Initialize and allocate memory for KINSOL

  kmem = FKINCreate(sunctx)
  if (.not. c_associated(kmem)) then
     print *, 'ERROR: kmem = NULL'
     stop 1
  end if  

  ! sunvec_u is used as a template

  ierr = FKINInit(kmem, c_funloc(func), sunvec_u)
  if (ierr /= 0) then
     print *, 'Error in FKINInit, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  ierr = FKINSetMaxSetupCalls(kmem, msbpre)
  if (ierr /= 0) then
     print *, 'Error in FKINSetMaxSetupCalls, ierr = ', ierr, '; halting'
     stop 1
  end if

  ftol = fnormtol
  ierr = FKINSetFuncNormTol(kmem, ftol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetFuncNormTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINSetScaledStepTol(kmem, scsteptol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetScaledStepTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINSetConstraints(kmem, constr)
  if (ierr /= 0) then
     print *, 'Error in FKINSetConstraints, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Create a sparse SPGMR linear solver

  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_u, prectype, maxl, sunctx)
  if (.not. associated(sunlinsol_LS)) then
     print *,'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! -------------------------
  ! Attach sparse linear solver

  ierr = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_J)
  if (ierr /= 0) then
     print *, 'Error in FKINSetLinearSolver, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set more optional SPGMR inputs

  ierr = FSUNLinSol_SPGMRSetMaxRestarts(kmem, maxlrst)
  if (ierr /= 0) then
     print *, 'Error in FSUNLinSol_SPGMRSetMaxRestarts, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FSUNLinSol_SPGMRSetPrecType(kmem, 1)
  if (ierr /= 0) then
     print *, 'Error in FSUNLinSol_SPGMRSetPrecType, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set preconditioner functions

  ierr = FKINSetPreconditioner(kmem, c_funloc(kpsetup), c_funloc(kpsolve))
  if (ierr /= 0) then
     print *, 'Error in FKINSetPreconditioner, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Call KINSol to solve problem
  !
  ! arguments: KINSol memory block
  !            Initial guess on input, solution on output
  !            Globalization strategy choice
  !            Scaling vector for the solution
  !            Scaling vector for the residual

  ierr = FKINSol(kmem, sunvec_u, KIN_NONE, sunvec_s, sunvec_s)
  if (ierr /= 0) then
     print *, 'Error in FKINSol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Print solution and solver statistics

  ! Get scaled norm of the system function
  ierr = FKINGetFuncNorm(kmem, fnorm)
  if (ierr /= 0) then
     print *, 'Error in FKINGetFuncNorm, ierr = ', ierr, '; halting'
     stop 1
  end if
  print *, " "
  print *, "Computed solution (||F|| = ", fnorm,"):"
  print *, " "
  call PrintOutput(u)
  call PrintFinalStats(kmem)

  ! clean up
  call FKINFree(kmem)
  ierr = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_J)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_s)
  ierr = FSUNContext_Free(sunctx)

end program main


! ----------------------------------------------------------------
! PrintOutput: prints solution at selected points
! ----------------------------------------------------------------
subroutine PrintOutput(uu)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use diag_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double), dimension(neq) :: uu

  !======= Internals ============

  do i = 1,neq,4
    print (i4, 4(1x, f10.6)), i, uu(i), uu(i+1), uu(i+2), uu(i+3)
  continue

  return

end subroutine PrintOutput


! ----------------------------------------------------------------
! PrintFinalStats
!
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(kmemo)

  !======= Inclusions ===========
  use iso_c_binding
  use fkinsol_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), parameter :: kmemo = kmem

  integer(c_long) :: nni(1), nli(1), nfe(1), npe(1), nps(1), ncfl(1)

  !======= Internals ============

  ! Main solver statistics

  ierr = FKINGetNumNonlinSolvIters(kmemo, nni)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumNonlinSolvIters, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumLinIters(kmemo, nli)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinIters, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumFuncEvals(kmemo, nfe)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = KINGetNumPrecEvals(kmemo, npe)
  if (ierr /= 0) then
    print *, 'Error in KINGetNumPrecEvals, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = KINGetNumPrecSolves(kmemo, nps)
  if (ierr /= 0) then
    print *, 'Error in KINGetNumPrecSolves, ierr = ', ierr, '; halting'
    stop 1
  end if

  ierr = KINGetNumLinConvFails(kmemo, nlcf)
  if (ierr /= 0) then
    print *, 'Error in KINGetNumLinConvFails, ierr = ', ierr, '; halting'
    stop 1
  end if

  print *, ' '
  print *, 'Final Statistics..'
  print *, ' '
  print '(2(A,i6))'    ,'nni      =', nni,      '    nli     =', nli
  print '(2(A,i6))'    ,'nfe      =', nfe,      '    npe     =', npe
  print '(2(A,i6))'    ,'nps      =', nps,      '    nlcf    =', nlcf

  return

end subroutine PrintFinalStats
