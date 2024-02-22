! ----------------------------------------------------------------
! Programmer(s): Daniel M. Margolis @ SMU
!                modified from previous ARKODE examples
! ----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2023, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ----------------------------------------------------------------
! Example problem:
!
! Diagonal ODE example. Stiff case, with BDF/SPGMR, diagonal
! preconditioner. Solved with preconditioning on left, then with
! preconditioning on right.
!
! ----------------------------------------------------------------

module DiagkryData
  ! --------------------------------------------------------------
  ! Description:
  !    Module containing problem-defining parameters.
  ! --------------------------------------------------------------

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  save

  ! SUNDIALS simulation context
  type(c_ptr) :: sunctx

  ! MPI domain decomposition information
  integer, target :: comm  ! communicator object
  integer :: myid          ! MPI process ID
  integer :: nprocs        ! total number of MPI processes

  ! Problem parameters
  integer(c_int),  parameter :: iGStype = 1
  integer(c_int),  parameter :: iPretype0 = 1
  integer(c_long), parameter :: nlocal = 10
  integer(c_long) :: neq
  integer(c_int) :: iPretype
  real(c_double) :: alpha

contains

  ! ----------------------------------------------------------------
  ! ODE RHS function f(t,y) (implicit).
  ! ----------------------------------------------------------------
  integer(c_int) function firhs(t, sunvec_y, sunvec_ydot, user_data) &
       result(retval) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t            ! current time
    type(N_Vector)        :: sunvec_y     ! solution N_Vector
    type(N_Vector)        :: sunvec_ydot  ! rhs N_Vector
    type(c_ptr),    value :: user_data    ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: y(:)
    real(c_double), pointer, dimension(nlocal) :: ydot(:)

    ! local data
    integer :: i

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    y(1:nlocal)    => FN_VGetArrayPointer(sunvec_y)
    ydot(1:nlocal) => FN_VGetArrayPointer(sunvec_ydot)

    ! Initialize ydot to zero
    ydot = 0.d0

    ! Fill ydot with rhs function
    do i = 1,nlocal
       ydot(i) = -alpha * (myid * nlocal + i) * y(i)
    end do

    retval = 0              ! Return with success
    return
  end function firhs
  ! ----------------------------------------------------------------

  ! ----------------------------------------------------------------
  ! Routine to solve preconditioner linear system
  ! This routine uses a diagonal preconditioner P = I - gamma*J,
  ! where J is a diagonal approximation to the true Jacobian, given by:
  ! J = diag(0, 0, 0, -4*alpha, ..., -N*alpha).
  ! The vector r is copied to z, and the inverse of P (restricted to the
  ! local vector segment) is applied to the vector z.
  ! ----------------------------------------------------------------
  integer(c_int) function Psolve(t, sunvec_y, sunvec_f, sunvec_r, sunvec_z, &
       gamma, delta, lr, user_data) result(retval) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t              ! current time
    type(N_Vector)        :: sunvec_y       ! solution N_Vector
    type(N_Vector)        :: sunvec_f       ! rhs N_Vector
    type(N_Vector)        :: sunvec_r       ! rhs N_Vector
    type(N_Vector)        :: sunvec_z       ! rhs N_Vector
    real(c_double), value :: gamma          ! current gamma value
    real(c_double), value :: delta          ! current delta value
    integer(c_int), value :: lr             ! left or right preconditioning
    type(c_ptr),    value :: user_data      ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer, dimension(nlocal) :: z(:)
    real(c_double), pointer, dimension(nlocal) :: r(:)

    ! local data
    integer(c_long) :: i, ibase, istart
    real(c_double)  :: psubi, pj

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    z(1:nlocal) => FN_VGetArrayPointer(sunvec_z)
    r(1:nlocal) => FN_VGetArrayPointer(sunvec_r)

    ! Copy r into z
    z = r

    ! Calculate Jacobian here
    ibase = myid * nlocal
    istart = max(1, 4 - ibase)
    do i = istart,nlocal
       pj = dble(ibase + i)
       psubi = 1.d0 + gamma * alpha * pj
       z(i) = z(i) / psubi
    end do

    retval = 0              ! Return with success
    return
  end function Psolve
  ! ----------------------------------------------------------------

end module DiagkryData
! ------------------------------------------------------------------


! ------------------------------------------------------------------
! Main driver program
! ------------------------------------------------------------------
program driver

  ! inclusions
  use, intrinsic :: iso_c_binding
  use fsundials_futils_mod       ! Fortran utilities
  use fcvode_mod                 ! Access CVode
  use fsundials_types_mod        ! sundials defined types
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod      ! Access generic N_Vector
  use fnvector_parallel_mod      ! Access parallel N_Vector
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use fsunlinsol_spgmr_mod       ! Fortran interface to spgmr SUNLinearSolver
  use fsundials_context_mod      ! Access sundials context
  use DiagkryData

  !======= Declarations =========
  implicit none

  ! Declarations
  ! general problem parameters
  integer, parameter :: Nt = 10                 ! total number of output times
  real(c_double), parameter :: T0 = 0.d0        ! initial time
  real(c_double), parameter :: Tf = 1.d0        ! final time
  real(c_double), parameter :: rtol = 1.d-5     ! relative and absolute tolerances
  real(c_double), parameter :: atol = 1.d-10

  ! solution vector and other local variables
  type(SUNLinearSolver), pointer :: sunls            ! sundials linear solver
  type(SUNMatrix), pointer :: sunmat_A               ! sundials matrix (empty)
  type(N_Vector), pointer :: sunvec_y                ! solution N_Vector
  real(c_double), pointer, dimension(nlocal) :: y(:) ! vector data
  type(c_ptr)     :: cvode_mem                       ! CVODE memory
  integer(c_long) :: N, Ntot
  integer(c_int) :: retval
  integer :: ierr
  logical :: outproc
  real(c_double) :: t(1), dTout, tout
  integer(c_long) :: nst(1)      ! number of time steps
  integer(c_long) :: nfe(1)      ! number of RHS evals
  integer(c_long) :: netf(1)     ! number of error test fails
  integer(c_long) :: nni(1)      ! number of nonlinear iters
  integer(c_long) :: ncfn(1)     ! number of nonlinear convergence fails
  integer(c_long) :: ncfl(1)     ! number of linear convergence fails
  integer(c_long) :: nli(1)      ! number of linear iters
  integer(c_long) :: npre(1)     ! number of preconditioner setups
  integer(c_long) :: npsol(1)    ! number of preconditioner solves
  integer(c_long) :: lenrw(1)    ! main solver real/int workspace size
  integer(c_long) :: leniw(1)
  integer(c_long) :: lenrwls(1)  ! linear solver real/int workspace size
  integer(c_long) :: leniwls(1)
  real(c_double)  :: avdim(1)    ! avg Krylov subspace dim (NLI/NNI)
  integer :: i, ioutput
  real(c_double) :: errmax, erri, gerrmax

  ! Initialize MPI variables
  comm = MPI_COMM_WORLD
  myid = 0
  nprocs = 0

  ! initialize MPI
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Init = ", ierr
     stop 1
  end if
  call MPI_Comm_size(comm, nprocs, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Comm_size = ", ierr
     call MPI_Abort(comm, 1, ierr)
  end if
  call MPI_Comm_rank(comm, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Comm_rank = ", ierr
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Set input arguments neq and alpha
  neq = nprocs * nlocal
  alpha = 10.0d0

  ! Create SUNDIALS simulation context, now that comm has been configured
  retval = FSUNContext_Create(comm, sunctx)
  if (retval /= 0) then
     print *, "Error: FSUNContext_Create returned ", retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Initial problem output
  outproc = (myid == 0)
  if (outproc) then
     write(6,*) "  "
     write(6,*) "Diagonal test problem:";
     write(6,'(A,i4)') "   neq = " , neq
     write(6,'(A,i4)') "   nlocal = " , nlocal
     write(6,'(A,i4)') "   nprocs = " , nprocs
     write(6,'(A,es9.2)') "   rtol = ", rtol
     write(6,'(A,es9.2)') "   atol = ", atol
     write(6,'(A,es9.2)') "   alpha = ", alpha
     write(6,*) "   ydot_i = -alpha*i * y_i (i = 1,...,neq)"
     write(6,*) "   Method is BDF/NEWTON/SPGMR"
     write(6,*) "   Diagonal preconditioner uses approximate Jacobian"
     write(6,*) "  "
  endif

  ! Create solution vector, point at its data, and set initial condition
  sunvec_y => FN_VNew_Parallel(comm, nlocal, neq, sunctx)
  y(1:nlocal) => FN_VGetArrayPointer(sunvec_y)
  y = 1.d0

  ! Create the CVode timestepper module
  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  if (.not. c_associated(cvode_mem)) then
     print *, "Error: FCVodeCreate returned NULL"
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FCVodeInit(cvode_mem, c_funloc(firhs), t0, sunvec_y)
  if (retval /= 0) then
     print *, "Error: FCVodeInit returned ", retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Tell CVODE to use a SPGMR linear solver.
  sunls => FSUNLinSol_SPGMR(sunvec_y, iPretype0, 0, sunctx)
  if (.not. associated(sunls)) then
     print *, 'ERROR: sunls = NULL'
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Attach the linear solver (with NULL SUNMatrix object)
  sunmat_A => null()
  retval = FCVodeSetLinearSolver(cvode_mem, sunls, sunmat_A)
  if (retval /= 0) then
     print *, 'Error in FCVodeSetLinearSolver, retval = ', retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FSUNLinSol_SPGMRSetGSType(sunls, iGStype)
  if (retval /= 0) then
     print *, 'Error in FSUNLinSol_SPGMRSetGSType, retval = ', retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (retval /= 0) then
     print *, "Error: FCVodeSStolerances returned ", retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FCVodeSetPreconditioner(cvode_mem, c_null_ptr, c_funloc(Psolve))
  if (retval /= 0) then
     print *, "Error: FCVodeSetPreconditioner returned ", retval
     call MPI_Abort(comm, 1, ierr)
  end if

  do iPretype = 1,2

     if (iPretype == 2) then

        y = 1.d0

        retval = FCVodeReInit(cvode_mem, t0, sunvec_y)
        if (retval /= 0) then
           print *, "Error in FCVodeReInit, retval = ", retval
           call MPI_Abort(comm, 1, ierr)
        end if

        retval = FSUNLinSol_SPGMRSetPrecType(sunls, iPretype)
        if (retval /= 0) then
           print *, "Error in FSUNLinSol_SPGMRSetPrecType, retval = ", retval
           call MPI_Abort(comm, 1, ierr)
        end if

        if (outproc) write(6,*) "   Preconditioning on right:"

     end if

     if (iPretype == 1 .and. outproc) write(6,*) "   Preconditioning on left:"

     ! Main time-stepping loop: calls CVode to perform the integration, then
     ! prints results.  Stops when the final time has been reached
     t(1) = T0
     dTout = 0.1d0
     tout = T0+dTout
     if (outproc) then
        write(6,*) "        t         steps     fe"
        write(6,*) "   --------------------------------"
     end if
     do ioutput=1,Nt

        ! Integrate to output time
        retval = FCVode(cvode_mem, tout, sunvec_y, t, CV_NORMAL)
        if (retval /= 0) then
           print *, "Error: FCVode returned ", retval
           call MPI_Abort(comm, 1, ierr)
        end if

        retval = FCVodeGetNumSteps(cvode_mem, nst)
        if (retval /= 0) then
           print *, "Error: FCVodeGetNumSteps returned ", retval
           call MPI_Abort(comm, 1, ierr)
        end if

        retval = FCVodeGetNumRhsEvals(cvode_mem, nfe)
        if (retval /= 0) then
           print *, "Error: FCVodeGetNumRhsEvals returned ", retval
           call MPI_Abort(comm, 1, ierr)
        end if

        ! print solution stats and update internal time
        if (outproc)   write(6,'(3x,f10.6,3(3x,i6))') t, nst, nfe
        tout = min(tout + dTout, Tf)

     end do
     if (outproc) then
        write(6,*) "   --------------------------------"
     end if

     ! Get max. absolute error in the local vector.
     errmax = 0.d0
     do i = 1,nlocal
        erri = y(i) - exp(-alpha * (myid * nlocal + i) * t(1))
        errmax = max(errmax, abs(erri))
     end do

     ! Get global max. error from MPI_Reduce call.
     call MPI_Reduce(errmax, gerrmax, 1, MPI_DOUBLE, MPI_MAX, &
          0, comm, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error in MPI_Reduce = ", ierr
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Print global max. error
     if (outproc) print '(a,es10.2)', "Max. absolute error is ", gerrmax

     ! Get final statistics
     retval = FCVodeGetNumSteps(cvode_mem, nst)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumSteps returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumRhsEvals(cvode_mem, nfe)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumRhsEvals returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumPrecEvals(cvode_mem, npre)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumPrecEvals returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumPrecSolves(cvode_mem, npsol)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumPrecSolves returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumNonlinSolvIters(cvode_mem, nni)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumNonlinSolvIters returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumLinIters(cvode_mem, nli)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumLinIters returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     avdim = dble(nli) / dble(nni)

     retval = FCVodeGetNumNonlinSolvConvFails(cvode_mem, ncfn)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumNonlinSolvConvFails returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumLinConvFails(cvode_mem, ncfl)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumLinSolvConvFails returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetNumErrTestFails(cvode_mem, netf)
     if (retval /= 0) then
        print *, "Error: FCVodeGetNumErrTestFails returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetWorkSpace(cvode_mem, lenrw, leniw)
     if (retval /= 0) then
        print *, "Error: FCVodeGetWorkSpace returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FCVodeGetLinWorkSpace(cvode_mem, lenrwls, leniwls)
     if (retval /= 0) then
        print *, "Error: FCVodeGetLinWorkSpace returned ", retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Print some final statistics
     if (outproc) then
        write(6,*) "  "
        write(6,*) "Final Solver Statistics:"
        write(6,'(A,i6)') "   Internal solver steps = ", nst
        write(6,'(A,i6)') "   Total RHS evals = ", nfe
        write(6,'(A,i6)') "   Total preconditioner setups = ", npre
        write(6,'(A,i6)') "   Total preconditioner solves = ", npsol
        write(6,'(A,i6)') "   Total nonlinear iterations  = ", nni
        write(6,'(A,i6)') "   Total linear iterations     = ", nli
        write(6,'(A,f8.4)') "   Average Krylov subspace dimension = ", avdim
        write(6,'(A,i6)') "   Total Convergence Failures - Nonlinear = ", ncfn
        write(6,'(A,i6)') "                              - Linear    = ", ncfl
        write(6,'(A,i6)') "   Total number of error test failures = ", netf
        write(6,'(A,2i6)') "   Main solver real/int workspace sizes = ", lenrw, leniw
        write(6,'(A,2i6)') "   Linear solver real/int workspace sizes = ", lenrwls, leniwls
        write(6,'(A)') "    "
        write(6,'(A)') "    "
        write(6,'(A)') "    "
     end if
  end do

  ! Clean up and return with successful completion
  call FCVodeFree(cvode_mem)       ! free integrator memory
  call FN_VDestroy(sunvec_y)          ! free vector memory
  call MPI_Barrier(comm, ierr)
  call MPI_Finalize(ierr)             ! Finalize MPI

end program driver
! ----------------------------------------------------------------
