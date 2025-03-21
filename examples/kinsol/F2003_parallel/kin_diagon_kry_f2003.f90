! ------------------------------------------------------------------
! Programmer(s): Allan G. Taylor, Alan C. Hindmarsh and Radu Serban @ LLNL
!                modified by Daniel M. Margolis @ SMU
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

module kinDiagonKry_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  integer(c_int64_t), parameter :: neq = 128

  integer(c_int)  :: ierr, retval, nprint
  integer(c_int64_t) :: i, nlocal
  real(c_double), pointer, dimension(neq) :: u(:), scale(:), constr(:)
  real(c_double)             :: p(neq)
  integer(c_int), parameter :: prectype = 2
  integer(c_int), parameter :: maxl = 10
  integer(c_int), parameter :: maxlrst = 2
  integer(c_long), parameter :: msbpre = 5
  real(c_double), parameter :: fnormtol = 1.0d-5
  real(c_double), parameter :: scsteptol = 1.0d-4

  ! MPI domain decomposition information
  integer, target :: comm  ! communicator object
  integer :: myid          ! MPI process ID
  integer :: nprocs        ! total number of MPI processes

contains

  ! ----------------------------------------------------------------
  ! init: Initializes variables u, scale, constr
  ! ----------------------------------------------------------------
  subroutine init(sunvec_u, sunvec_s, sunvec_c)

    !======= Inclusions ===========
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_u  ! solution N_Vector
    type(N_Vector)       :: sunvec_s  ! scaling N_Vector
    type(N_Vector)       :: sunvec_c  ! constraint N_Vector

    ! local variables
    integer(c_int64_t) :: ii

    u(1:nlocal) => FN_VGetArrayPointer(sunvec_u)
    scale(1:nlocal) => FN_VGetArrayPointer(sunvec_s)
    constr(1:nlocal) => FN_VGetArrayPointer(sunvec_c)

    ! -------------------------
    ! Set initial guess, and disable scaling

    do i = 1, nlocal
      ii = i + myid*nlocal
      u(i) = 2.0d0*dble(ii)
    end do
    scale = 1.0d0
    constr = 0.0d0

  end subroutine init
  ! ----------------------------------------------------------------

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
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: uu(:), ff(:)

    ! local variables
    integer(c_int64_t) :: ii

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    uu(1:nlocal) => FN_VGetArrayPointer(sunvec_u)
    ff(1:nlocal) => FN_VGetArrayPointer(sunvec_f)

    ! loop over domain, computing our system f(u) = 0
    do i = 1, nlocal
      ! set local variables
      ii = i + myid*nlocal

      ! applying the constraint f(u) = u(i)^2 - i^2
      ff(i) = uu(i)*uu(i) - dble(ii*ii)
    end do

    ! return success
    ierr = 0
    return

  end function func
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! kpsetup: The KINSOL Preconditioner setup function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function kpsetup(sunvec_u, sunvec_s, sunvec_f, &
                                  sunvec_fs, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_s  ! scaling N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(N_Vector)     :: sunvec_fs ! LHS scaling N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: udata(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    udata(1:nlocal) => FN_VGetArrayPointer(sunvec_u)

    ! loop over domain
    do i = 1, nlocal

      ! setup preconditioner
      p(i) = 0.5d0/(udata(i) + 5.0d0)
    end do

    ! return success
    ierr = 0
    return

  end function kpsetup
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! kpsolve: The KINSOL Preconditioner solve function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function kpsolve(sunvec_u, sunvec_s, sunvec_f, &
                                  sunvec_fs, sunvec_v, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use fsundials_core_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_u  ! solution N_Vector
    type(N_Vector)     :: sunvec_s  ! scaling N_Vector
    type(N_Vector)     :: sunvec_f  ! LHS N_Vector
    type(N_Vector)     :: sunvec_fs ! LHS scaling N_Vector
    type(N_Vector)     :: sunvec_v  ! LHS scaling N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: v(:)

    !======= Internals ============
    ! get data arrays from SUNDIALS vectors
    v(1:nlocal) => FN_VGetArrayPointer(sunvec_v)

    ! loop over domain
    do i = 1, nlocal

      ! preconditioner solver
      v(i) = v(i)*p(i)
    end do

    ! return success
    ierr = 0
    return

  end function kpsolve
  ! ----------------------------------------------------------------

end module kinDiagonKry_mod
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use fsundials_core_mod
  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_parallel_mod      ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to SPGMR SUNLinearSolver
  use kinDiagonKry_mod           ! problem-defining functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: ftol

  type(c_ptr)                    :: sunctx        ! sundials context
  type(N_Vector), pointer :: sunvec_u      ! sundials vectors
  type(N_Vector), pointer :: sunvec_s      ! sundials vectors
  type(N_Vector), pointer :: sunvec_c      ! sundials vectors
  type(SUNMatrix), pointer :: sunmat_J      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver

  type(c_ptr) :: kmem ! KINSOL memory
  logical :: outproc

  !======= Internals ============

  ! -------------------------
  ! Initialize MPI variables
  comm = MPI_COMM_WORLD
  myid = 0
  nprocs = 0

  ! initialize MPI
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Init = ", ierr
    stop 1
  end if
  call MPI_Comm_size(comm, nprocs, ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Comm_size = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if
  if (popcnt(nprocs) /= 1 .or. nprocs > neq) then
    write (0, *) "Error nprocs must equal a power of 2^n <= neq for functionality."
    call MPI_Abort(comm, 1, ierr)
  end if
  call MPI_Comm_rank(comm, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
    write (0, *) "Error in MPI_Comm_rank = ", ierr
    call MPI_Abort(comm, 1, ierr)
  end if

  outproc = (myid == 0)

  ! Print problem description

  if (outproc) then
    print *, " "
    print *, "Example program kinDiagon_kry_f2003:"
    print *, "   This FKINSOL example solves a 128 eqn diagonal algebraic system."
    print *, " Its purpose is to demonstrate the use of the Fortran interface in"
    print *, " a parallel environment."
    print *, " "
    print *, "Solution method: KIN_none"
    print '(a,i3)', "Problem size: neq = ", neq
    print '(a,i3)', "Number of procs: nprocs = ", nprocs
  end if

  ! -------------------------
  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
  if (retval /= 0) then
    print *, 'ERROR in FSUNContext_Create'
    stop 1
  end if

  ! -------------------------
  ! Create vectors for solution and scales
  nlocal = neq/nprocs

  sunvec_u => FN_VNew_Parallel(comm, nlocal, neq, sunctx)
  sunvec_s => FN_VNew_Parallel(comm, nlocal, neq, sunctx)
  sunvec_c => FN_VNew_Parallel(comm, nlocal, neq, sunctx)

  call init(sunvec_u, sunvec_s, sunvec_c)

  ! -------------------------
  ! Initialize and allocate memory for KINSOL

  kmem = FKINCreate(sunctx)
  if (.not. c_associated(kmem)) then
    print *, 'ERROR: kmem = NULL'
    stop 1
  end if

  ! sunvec_u is used as a template

  retval = FKINInit(kmem, c_funloc(func), sunvec_u)
  if (retval /= 0) then
    print *, 'Error in FKINInit, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  retval = FKINSetMaxSetupCalls(kmem, msbpre)
  if (retval /= 0) then
    print *, 'Error in FKINSetMaxSetupCalls, retval = ', retval, '; halting'
    stop 1
  end if

  ftol = fnormtol
  retval = FKINSetFuncNormTol(kmem, ftol)
  if (retval /= 0) then
    print *, 'Error in FKINSetFuncNormTol, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINSetScaledStepTol(kmem, scsteptol)
  if (retval /= 0) then
    print *, 'Error in FKINSetScaledStepTol, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINSetConstraints(kmem, sunvec_c)
  if (retval /= 0) then
    print *, 'Error in FKINSetConstraints, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Create a SPGMR linear solver

  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_u, prectype, maxl, sunctx)
  if (.not. associated(sunlinsol_LS)) then
    print *, 'ERROR: sunlinsol = NULL'
    stop 1
  end if

  ! -------------------------
  ! Attach linear solver

  sunmat_J => null()

  retval = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_J)
  if (retval /= 0) then
    print *, 'Error in FKINSetLinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set more optional SPGMR inputs

  retval = FSUNLinSol_SPGMRSetMaxRestarts(sunlinsol_LS, maxlrst)
  if (retval /= 0) then
    print *, 'Error in FSUNLinSol_SPGMRSetMaxRestarts, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Set preconditioner functions

  retval = FKINSetPreconditioner(kmem, c_funloc(kpsetup), c_funloc(kpsolve))
  if (retval /= 0) then
    print *, 'Error in FKINSetPreconditioner, retval = ', retval, '; halting'
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

  retval = FKINSol(kmem, sunvec_u, KIN_NONE, sunvec_s, sunvec_s)
  if (retval /= 0) then
    print *, 'Error in FKINSol, retval = ', retval, '; halting'
    stop 1
  end if

  ! -------------------------
  ! Print solution and solver statistics

  if (outproc) then
    print *, " "
  end if
  do nprint = 0, nprocs - 1
    if (nprint == myid) then
      call PrintOutput(u)
    end if
    call MPI_Barrier(comm, ierr)
  end do
  call MPI_Barrier(comm, ierr)
  call PrintFinalStats(kmem, outproc)

  ! clean up
  call FKINFree(kmem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_s)
  call FN_VDestroy(sunvec_c)
  retval = FSUNContext_Free(sunctx)
  call MPI_Barrier(comm, ierr)
  call MPI_Finalize(ierr)             ! Finalize MPI

end program main
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintOutput: prints solution at selected points
! ----------------------------------------------------------------
subroutine PrintOutput(uu)

  !======= Inclusions ===========
  use kinDiagonKry_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double), dimension(neq) :: uu
  integer(c_int64_t) :: ii

  !======= Internals ============

  do i = 1, nlocal, 4
    ii = i + nlocal*myid
    print '(i4, 4(1x, f10.6))', ii, uu(i), uu(i + 1), uu(i + 2), uu(i + 3)
  end do

  return

end subroutine PrintOutput
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! PrintFinalStats
!
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(kmemo, outproc)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fkinsol_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: kmemo
  logical, intent(in) :: outproc

  integer(c_int) :: retval
  integer(c_long) :: nni(1), nli(1), nfe(1), npe(1), nps(1), ncfl(1)

  !======= Internals ============

  ! Main solver statistics

  retval = FKINGetNumNonlinSolvIters(kmemo, nni)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumNonlinSolvIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumLinIters(kmemo, nli)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumLinIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumFuncEvals(kmemo, nfe)
  if (retval /= 0) then
    print *, 'Error in FKINGetNumFuncEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumPrecEvals(kmemo, npe)
  if (retval /= 0) then
    print *, 'Error in KINGetNumPrecEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumPrecSolves(kmemo, nps)
  if (retval /= 0) then
    print *, 'Error in KINGetNumPrecSolves, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FKINGetNumLinConvFails(kmemo, ncfl)
  if (retval /= 0) then
    print *, 'Error in KINGetNumLinConvFails, retval = ', retval, '; halting'
    stop 1
  end if

  if (outproc) then
    print *, ' '
    print *, 'Final Statistics..'
    print *, ' '
    print '(2(A,i6))', 'nni      =', nni, '    nli     =', nli
    print '(2(A,i6))', 'nfe      =', nfe, '    npe     =', npe
    print '(2(A,i6))', 'nps      =', nps, '    nlcf    =', ncfl
  end if

  return

end subroutine PrintFinalStats
! ----------------------------------------------------------------
