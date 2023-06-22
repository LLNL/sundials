!-----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
!                Daniel M. Margolis @ SMU
!-----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2023, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
!-----------------------------------------------------------------
! Example problem:
!
! The following test simulates a simple anisotropic 2D heat
! equation,
!    u_t = kx*u_xx + ky*u_yy + h,
! for t in [0, 10], (x,y) in [0, 1]^2, with initial conditions
!    u(0,x,y) =  0,
! stationary boundary conditions, i.e.
!    u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
! and a heat source of the form
!    h(x,y) = sin(pi*x)*sin(2*pi*y).
!
! Under this setup, the problem has an analytical solution:
!    u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
!    a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).
!
! The spatial derivatives are computed using second-order
! centered differences, with the data distributed over nx*ny
! points on a uniform spatial grid.
!
! The spatial grid parameters nx and ny, the parameters kx and ky,
! as well as the desired relative and absolute solver tolerances,
! are provided in the input file input_heat2D.txt.
!
! This program solves the problem with a DIRK method.  This
! employs a Newton iteration with the PCG iterative linear solver,
! which itself uses a Jacobi preconditioner.  The example uses the
! built-in finite-difference Jacobian-vector product routine, but
! supplies both the RHS and preconditioner setup/solve functions.
!
! 20 outputs are printed at equal intervals, and run statistics
! are printed at the end.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
module HeatUserData
  !---------------------------------------------------------------
  ! Description:
  !    Module containing problem-defining parameters, as well as
  !    data buffers for MPI exchanges with neighboring processes.
  !    Also contains routines to:
  !      (a) initialize the module
  !      (b) perform exchanges
  !      (c) free module data.
  !---------------------------------------------------------------
   use, intrinsic :: iso_c_binding
  implicit none
  save

  integer(c_long)  :: nx                                ! global number of x grid points
  integer(c_long)  :: ny                                ! global number of y grid points
  integer(c_long)  :: is                                ! global x indices of this subdomain
  integer(c_long)  :: ie
  integer(c_long)  :: js                                ! global y indices of this subdomain
  integer(c_long)  :: je
  integer(c_long)  :: nxl                               ! local number of x grid points
  integer(c_long)  :: nyl                               ! local number of y grid points
  double precision :: dx                                ! x-directional mesh spacing
  double precision :: dy                                ! y-directional mesh spacing
  double precision :: kx                                ! x-directional diffusion coefficient
  double precision :: ky                                ! y-directional diffusion coefficient
  double precision, dimension(:,:), allocatable :: h    ! heat source vector
  double precision, dimension(:,:), allocatable :: d    ! inverse of Jacobian diagonal
  integer, target  :: comm                              ! communicator object
  integer(c_int)   :: myid                              ! MPI process ID
  integer(c_int)   :: nprocs                            ! total number of MPI processes
!   type(c_ptr) :: sunctx                                 ! SUNDIALS simulation context
!   type(c_ptr) :: logger                                 ! SUNDIALS logger (optional, must change CMake)
  logical :: HaveNbor(2,2)                              ! flags denoting neighbor on boundary
  double precision, dimension(:), allocatable :: Erecv  ! receive buffers for neighbor exchange
  double precision, dimension(:), allocatable :: Wrecv
  double precision, dimension(:), allocatable :: Nrecv
  double precision, dimension(:), allocatable :: Srecv
  double precision, dimension(:), allocatable :: Esend  ! send buffers for neighbor exchange
  double precision, dimension(:), allocatable :: Wsend
  double precision, dimension(:), allocatable :: Nsend
  double precision, dimension(:), allocatable :: Ssend
  integer(c_long) :: nst(1), nst_att(1), nfe(1), nfi(1), nls(1), nli(1), njve(1), npe(1), nps(1), &
                           nlcf(1), niters(1), nlscf(1), netf(1)

contains

  !---------------------------------------------------------------
  ! Initialize memory allocated within Userdata (set to defaults) -- keep this around (with "type" adjustments)
  !---------------------------------------------------------------
  subroutine InitUserData()
    implicit none
    include "mpif.h"
    nx = 0
    ny = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    nxl = 0
    nyl = 0
    dx = 0.d0
    dy = 0.d0
    kx = 0.d0
    ky = 0.d0
    if (allocated(h))  deallocate(h)
    if (allocated(d))  deallocate(d)
    myid = 0
    nprocs = 0
    HaveNbor = .false.
    if (allocated(Erecv))  deallocate(Erecv)
    if (allocated(Wrecv))  deallocate(Wrecv)
    if (allocated(Nrecv))  deallocate(Nrecv)
    if (allocated(Srecv))  deallocate(Srecv)
    if (allocated(Esend))  deallocate(Esend)
    if (allocated(Wsend))  deallocate(Wsend)
    if (allocated(Nsend))  deallocate(Nsend)
    if (allocated(Ssend))  deallocate(Ssend)
  end subroutine InitUserData
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Set up parallel decomposition
  !---------------------------------------------------------------
  subroutine SetupDecomp(ierr)
    ! declarations
    implicit none
    include "mpif.h"
    integer(c_int), intent(inout) :: ierr
    integer(c_int) :: dims(2), periods(2), coords(2)

    ! internals

    ! check that this has not been called before
    if (allocated(h) .or. allocated(d)) then
       write(0,*) "SetupDecomp warning: parallel decomposition already set up"
       ierr = 1
       return
    end if

    ! get suggested parallel decomposition
    dims = (/0, 0/)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Comm_size = " , ierr
       return
    end if
    call MPI_Dims_create(nprocs, 2, dims, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Dims_create = " , ierr
       return
    end if

    ! set up 2D Cartesian communicator
    periods = (/0, 0/)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Cart_create = " , ierr
       return
    end if
    call MPI_Comm_rank(comm, myid, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Comm_rank = " , ierr
       return
    end if

    ! determine local extents
    call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Cart_get = " , ierr
       return
    end if
    is = nx*coords(1)/dims(1) + 1
    ie = nx*(coords(1)+1)/dims(1)
    js = ny*coords(2)/dims(2) + 1
    je = ny*(coords(2)+1)/dims(2)
    nxl = ie-is+1
    nyl = je-js+1

    ! determine if I have neighbors, and allocate exchange buffers
    HaveNbor(1,1) = (is /= 1)
    HaveNbor(1,2) = (ie /= nx)
    HaveNbor(2,1) = (js /= 1)
    HaveNbor(2,2) = (je /= ny)
    if (HaveNbor(1,1)) then
       allocate(Wrecv(nyl))
       allocate(Wsend(nyl))
    end if
    if (HaveNbor(1,2)) then
       allocate(Erecv(nyl))
       allocate(Esend(nyl))
    end if
    if (HaveNbor(2,1)) then
       allocate(Srecv(nxl))
       allocate(Ssend(nxl))
    end if
    if (HaveNbor(2,2)) then
       allocate(Nrecv(nxl))
       allocate(Nsend(nxl))
    end if

    ! allocate temporary vectors
    allocate(h(nxl,nyl))    ! Create vector for heat source
    allocate(d(nxl,nyl))    ! Create vector for Jacobian diagonal

    ierr = 0     ! return with success flag
    return
  end subroutine SetupDecomp
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Perform neighbor exchange -- keep this around, but you'll need to adjust "y" to have the correct type, etc.
  !---------------------------------------------------------------
  subroutine Exchange(sunvec_y, ierr)
      ! declarations
      use fsundials_nvector_mod
      implicit none
      include "mpif.h"
      type(N_Vector), intent(inout) :: sunvec_y
      integer(c_int), intent(out)   :: ierr
      real(c_double), pointer, dimension(:,:) :: y
      integer(c_int)  :: reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
      integer(c_int)  :: stat(MPI_STATUS_SIZE)
      integer(c_long) :: i
      integer(c_int)  :: ipW, ipE, ipS, ipN
      integer(c_int)  :: coords(2), dims(2), periods(2), nbcoords(2)

      ! internals

      y(1:nxl, 1:nyl) => FN_VGetArrayPointer(sunvec_y)

      ! MPI neighborhood information
      call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
      if (ierr /= MPI_SUCCESS) then
         write(0,*) "Error in MPI_Cart_get = ", ierr
         return
      end if
      if (HaveNbor(1,1)) then
         nbcoords = (/ coords(1)-1, coords(2) /)
         call MPI_Cart_rank(comm, nbcoords, ipW, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Cart_rank = ", ierr
            return
         end if
      end if
      if (HaveNbor(1,2)) then
         nbcoords = (/ coords(1)+1, coords(2) /)
         call MPI_Cart_rank(comm, nbcoords, ipE, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Cart_rank = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,1)) then
         nbcoords = (/ coords(1), coords(2)-1 /)
         call MPI_Cart_rank(comm, nbcoords, ipS, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Cart_rank = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,2)) then
         nbcoords = (/ coords(1), coords(2)+1 /)
         call MPI_Cart_rank(comm, nbcoords, ipN, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Cart_rank = ", ierr
            return
         end if
      end if

      ! open Irecv buffers
      if (HaveNbor(1,1)) then
         call MPI_Irecv(Wrecv, nyl, MPI_DOUBLE_PRECISION, ipW, &
                        MPI_ANY_TAG, comm, reqRW, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Irecv = ", ierr
            return
         end if
      end if
      if (HaveNbor(1,2)) then
         call MPI_Irecv(Erecv, nyl, MPI_DOUBLE_PRECISION, ipE, &
                        MPI_ANY_TAG, comm, reqRE, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Irecv = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,1)) then
         call MPI_Irecv(Srecv, nxl, MPI_DOUBLE_PRECISION, ipS, &
                        MPI_ANY_TAG, comm, reqRS, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Irecv = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,2)) then
         call MPI_Irecv(Nrecv, nxl, MPI_DOUBLE_PRECISION, ipN, &
                        MPI_ANY_TAG, comm, reqRN, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Irecv = ", ierr
            return
         end if
      end if

      ! send data
      if (HaveNbor(1,1)) then
         do i=1,nyl
            Wsend(i) = y(1,i)
         end do
         call MPI_Isend(Wsend, nyl, MPI_DOUBLE_PRECISION, ipW, 0, &
                        comm, reqSW, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Isend = ", ierr
            return
         end if
      end if
      if (HaveNbor(1,2)) then
         do i=1,nyl
            Esend(i) = y(nxl,i)
         end do
         call MPI_Isend(Esend, nyl, MPI_DOUBLE_PRECISION, ipE, 1, &
                        comm, reqSE, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Isend = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,1)) then
         do i=1,nxl
            Ssend(i) = y(i,1)
         end do
         call MPI_Isend(Ssend, nxl, MPI_DOUBLE_PRECISION, ipS, 2, &
                        comm, reqSS, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Isend = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,2)) then
         do i=1,nxl
            Nsend(i) = y(i,nyl)
         end do
         call MPI_Isend(Nsend, nxl, MPI_DOUBLE_PRECISION, ipN, 3, &
                        comm, reqSN, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Isend = ", ierr
            return
         end if
      end if

      ! wait for messages to finish
      if (HaveNbor(1,1)) then
         call MPI_Wait(reqRW, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
         call MPI_Wait(reqSW, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
      end if
      if (HaveNbor(1,2)) then
         call MPI_Wait(reqRE, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
         call MPI_Wait(reqSE, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,1)) then
         call MPI_Wait(reqRS, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
         call MPI_Wait(reqSS, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
      end if
      if (HaveNbor(2,2)) then
         call MPI_Wait(reqRN, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
         call MPI_Wait(reqSN, stat, ierr)
         if (ierr /= MPI_SUCCESS) then
            write(0,*) "Error in MPI_Wait = ", ierr
            return
         end if
      end if

      ierr = 0     ! return with success flag
      return
   end subroutine Exchange
   !---------------------------------------------------------------


   !---------------------------------------------------------------
   ! Free memory allocated within Userdata
   !---------------------------------------------------------------
   subroutine FreeUserData(ierr)
      implicit none
      integer(c_int), intent(out) :: ierr
      if (allocated(h))      deallocate(h)
      if (allocated(d))      deallocate(d)
      if (allocated(Wrecv))  deallocate(Wrecv)
      if (allocated(Wsend))  deallocate(Wsend)
      if (allocated(Erecv))  deallocate(Erecv)
      if (allocated(Esend))  deallocate(Esend)
      if (allocated(Srecv))  deallocate(Srecv)
      if (allocated(Ssend))  deallocate(Ssend)
      if (allocated(Nrecv))  deallocate(Nrecv)
      if (allocated(Nsend))  deallocate(Nsend)
      ierr = 0     ! return with success flag
      return
   end subroutine FreeUserData
   !---------------------------------------------------------------


   !---------------------------------------------------------------
   ! RMS norm function for parallel array data
   !---------------------------------------------------------------
   subroutine PRMS(sunvec_y,yrms,ierr)
      ! declarations
      use fsundials_nvector_mod
      implicit none
      include "mpif.h"
      type(N_Vector), intent(inout) :: sunvec_y
      integer(c_int), intent(out)   :: ierr
      real(c_double), pointer     :: y(:,:)
      double precision, intent(out) :: yrms
      double precision :: lsum, gsum

      ! internals
      y(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_y)
      lsum = sum(y**2)
      call MPI_Allreduce(lsum, gsum, 1, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, comm, ierr)
      if (ierr /= MPI_SUCCESS) then
         write(0,*) "Error in MPI_Allreduce = ", ierr
         call MPI_Finalize(ierr)
      end if
      yrms = sqrt(gsum/nx/ny)

      ierr = 0     ! return with success flag
      return
   end subroutine PRMS
   !---------------------------------------------------------------

   end module HeatUserData
   !-----------------------------------------------------------------

   module Implicit_and_Prec_Fn

      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      implicit none
      save

   contains

      !-----------------------------------------------------------------
      ! Function called by the solver
      !-----------------------------------------------------------------


      integer(c_int) function farkifun(t, sunvec_y, sunvec_ydot, user_data) &
               result(ierr) bind(C)
      !-----------------------------------------------------------------
      ! f routine to compute the ODE RHS function f(t,y)
      !-----------------------------------------------------------------
      ! declarations
      use HeatUserData

      implicit none
      include "mpif.h"

      real(c_double) :: t
      type(N_Vector)   :: sunvec_y
      type(N_Vector)   :: sunvec_ydot
      type(c_ptr)      :: user_data
      double precision :: c1, c2, c3
      integer(c_long)  :: ii, jj
      real(c_double), pointer, dimension(nxl,nyl) :: y(:,:)
      real(c_double), pointer, dimension(nxl,nyl) :: ydot(:,:)

      ! internals

      ! Initialize ydot to zero
      y(1:nxl,1:nyl)    => FN_VGetArrayPointer(sunvec_y)
      ydot(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_ydot)
      ydot = 0.d0

      ! Exchange boundary data with neighbors
      call Exchange(sunvec_y, ierr)
      if (ierr /= MPI_SUCCESS) then
         write(0,*) "Error in Exchange = " , ierr
         return
      end if

      ! iterate over subdomain interior, computing approximation to RHS
      c1 = kx/dx/dx
      c2 = ky/dy/dy
      c3 = -2.d0*(c1 + c2)
      do jj=2,nyl-1
         do ii=2,nxl-1
            ydot(ii,jj) = c1*(y(ii-1,jj)+y(ii+1,jj)) + c2*(y(ii,jj-1)+y(ii,jj+1)) + c3*y(ii,jj)
         end do
      end do

      ! iterate over subdomain boundaries (if not at overall domain boundary)
      if (HaveNbor(1,1)) then    ! West face
         ii=1
         do jj=2,nyl-1
            ydot(ii,jj) = c1*(Wrecv(jj)+y(ii+1,jj)) + c2*(y(ii,jj-1)+y(ii,jj+1)) + c3*y(ii,jj)
         end do
      end if
      if (HaveNbor(1,2)) then    ! East face
         ii=nxl
         do jj=2,nyl-1
            ydot(ii,jj) = c1*(y(ii-1,jj)+Erecv(jj)) + c2*(y(ii,jj-1)+y(ii,jj+1)) + c3*y(ii,jj)
         end do
      end if
      if (HaveNbor(2,1)) then    ! South face
         jj=1
         do ii=2,nxl-1
            ydot(ii,jj) = c1*(y(ii-1,jj)+y(ii+1,jj)) + c2*(Srecv(ii)+y(ii,jj+1)) + c3*y(ii,jj)
         end do
      end if
      if (HaveNbor(2,2)) then    ! West face
         jj=nyl
         do ii=2,nxl-1
            ydot(ii,jj) = c1*(y(ii-1,jj)+y(ii+1,jj)) + c2*(y(ii,jj-1)+Nrecv(ii)) + c3*y(ii,jj)
         end do
      end if
      if (HaveNbor(1,1) .and. HaveNbor(2,1)) then  ! South-West corner
         ii=1
         jj=1
         ydot(ii,jj) = c1*(Wrecv(jj)+y(ii+1,jj)) + c2*(Srecv(ii)+y(ii,jj+1)) + c3*y(ii,jj)
      end if
      if (HaveNbor(1,1) .and. HaveNbor(2,2)) then  ! North-West corner
         ii=1
         jj=nyl
         ydot(ii,jj) = c1*(Wrecv(jj)+y(ii+1,jj)) + c2*(y(ii,jj-1)+Nrecv(ii)) + c3*y(ii,jj)
      end if
      if (HaveNbor(1,2) .and. HaveNbor(2,1)) then  ! South-East corner
         ii=nxl
         jj=1
         ydot(ii,jj) = c1*(y(ii-1,jj)+Erecv(jj)) + c2*(Srecv(ii)+y(ii,jj+1)) + c3*y(ii,jj)
      end if
      if (HaveNbor(1,2) .and. HaveNbor(2,2)) then  ! North-East corner
         ii=nxl
         jj=nyl
         ydot(ii,jj) = c1*(y(ii-1,jj)+Erecv(jj)) + c2*(y(ii,jj-1)+Nrecv(ii)) + c3*y(ii,jj)
      end if

      ydot = ydot + h         ! add in heat source

      ierr = 0                ! Return with success
      return

      end function farkifun
      !-----------------------------------------------------------------


      integer(c_int) function farkpset(t, sunvec_y, sunvec_fy, jok, jcur, &
      gamma, user_data) result(ierr) bind(C)
      !-----------------------------------------------------------------
      ! Preconditioner setup routine (fills inverse of Jacobian diagonal)
      !-----------------------------------------------------------------
      ! declarations
      use HeatUserData
      implicit none
      type(N_Vector)   :: sunvec_y
      type(N_Vector)   :: sunvec_fy
      real(c_double) :: t, gamma
      double precision :: c
      integer(c_int), value :: jok
      integer(c_int)   :: jcur
      type(c_ptr)      :: user_data

      ! internals
      c = 1.d0 + gamma*2.d0*(kx/dx/dx + ky/dy/dy)

      ! set all entries of d to the inverse of the diagonal values in interior
      ! (since boundary RHS is 0, set boundary diagonals to the same)
      d = 1.d0/c

      jcur = 1 ! Update jcur flag
      ierr = 0 ! Return with success
      return
      end function farkpset


      integer(c_int) function farkpsol(t, sunvec_y, sunvec_fy, sunvec_r, sunvec_z, &
               gamma, delta, lr, user_data) result(ierr) bind(C)
      !-----------------------------------------------------------------
      ! Preconditioner solve routine
      !-----------------------------------------------------------------
      ! declarations
      use HeatUserData
      implicit none

      type(N_Vector)   :: sunvec_y, sunvec_fy, sunvec_r, sunvec_z
      real(c_double)   :: t, gamma, delta
      integer(c_int)   :: lr
      type(c_ptr)      :: user_data
      real(c_double), pointer, dimension(nxl, nyl) :: r(:,:), z(:,:)

      ! internals
      r(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_r)
      z(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_z)

      z = r*d      ! perform Jacobi iteration (whole array operation)
      ierr = 0     ! Return with success
      return
   end function farkpsol
   !-----------------------------------------------------------------

   subroutine GetFinalStats(arkmem,flag,retval)

      use HeatUserData
      use farkode_mod                   ! Access ARKode
      use farkode_arkstep_mod           ! Access ARKStep

      implicit none
      include "mpif.h"

      type(c_ptr),          intent(in) :: arkmem
      integer(c_int),       intent(in) :: flag
      integer(c_int),    intent(inout) :: retval

      retval = FARKStepGetNumSteps(arkmem, nst)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumSteps = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumStepAttempts(arkmem, nst_att)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumStepAttempts = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumRhsEvals(arkmem, nfe, nfi)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumRhsEvals = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumLinSolvSetups(arkmem, nls)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumLinSolvSetups = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumLinIters(arkmem, nli)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumLinIters = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumJtimesEvals(arkmem, njve)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumJtimesEvals = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumPrecEvals(arkmem, npe)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumPrecEvals = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumPrecSolves(arkmem, nps)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumPrecSolves = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumLinConvFails(arkmem, nlcf)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumLinConvFails = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumNonlinSolvIters(arkmem, niters)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumNonlinSolvIters = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumNonlinSolvConvFails(arkmem, nlscf)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumNonlinSolvConvFails = ", retval
         call MPI_Finalize(flag)
      end if
      retval = FARKStepGetNumErrTestFails(arkmem, netf)
      if (retval /= 0) then
         write(0,*) "Error in FARKStepGetNumErrTestFails = ", retval
         call MPI_Finalize(flag)
      end if

   end subroutine GetFinalStats

end module Implicit_and_Prec_Fn

module Output_Fns

      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod         ! Access generic N_Vector

      implicit none
      save

      integer(c_long) :: i, j

   contains

      subroutine InitOutput_Disk(sunvec_y, stream, nxl, nyl)

         implicit none

         type(N_Vector), intent(inout) :: sunvec_y
         integer(c_int), intent(in)    :: stream
         integer(c_long), intent(in)   :: nxl, nyl

         real(c_double), pointer, dimension(:,:) :: y
         y(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_y)

         ! Output initial condition to disk
         do j=1,nyl
            do i=1,nxl
               write(stream,'(es25.16)',advance='no') y(i,j)
            end do
         end do
         write(stream,*) "  "

      end subroutine InitOutput_Disk

      subroutine OutputResults_Disk(sunvec_y, stream, nxl, nyl)

         implicit none

         type(N_Vector), intent(inout) :: sunvec_y
         integer(c_int), intent(in)    :: stream
         integer(c_long), intent(in)   :: nxl, nyl

         real(c_double), pointer, dimension(:,:) :: y
         y(1:nxl,1:nyl) => FN_VGetArrayPointer(sunvec_y)

         ! output results to disk
         do j=1,nyl
            do i=1,nxl
               write(stream,'(es25.16)',advance='no') y(i,j)
            end do
         end do
         write(stream,*) "  "

      end subroutine OutputResults_Disk

end module Output_Fns

!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
program driver

  ! inclusions
  use HeatUserData
  use Implicit_and_Prec_Fn
  use Output_Fns
  use fsundials_types_mod           ! sundials defined types
  use fsundials_context_mod         ! Access sundials context
  use farkode_mod                   ! Access ARKode
  use farkode_arkstep_mod           ! Access ARKStep
  use fsundials_nvector_mod         ! Access generic N_Vector
  use fnvector_parallel_mod         ! Access parallel N_Vector
  use fsundials_matrix_mod          ! Access generic SUNMatrix
  use fsundials_linearsolver_mod    ! Access generic SUNLinearSolver
  use fsunlinsol_pcg_mod            ! Access PCG SUNLinearSolver

  implicit none
  include "mpif.h"

  ! Declarations
  ! general problem parameters
  double precision, parameter :: pi = 3.1415926535897932d0
  integer(c_int),   parameter :: Nt = 20           ! total number of output times
  integer(c_long),  parameter :: nx_ = 60          ! spatial mesh size
  integer(c_long),  parameter :: ny_ = 120
  integer(c_int),   parameter :: PCGpretype = 1    ! enable preconditioner
  integer(c_int),   parameter :: PCGmaxl = 20      ! max num. PCG iterations
  double precision, parameter :: T0 = 0.d0         ! initial time
  double precision, parameter :: Tf = 0.3d0        ! final time
  double precision, parameter :: rtol = 1.d-5      ! relative and absolute tolerances
  double precision, parameter :: atol = 1.d-10
  double precision, parameter :: kx_ = 0.5d0       ! heat conductivity coefficients
  double precision, parameter :: ky_ = 0.75d0
  double precision, parameter :: nlscoef = 1.d-7   ! nonlinear solver tolerance factor

  ! solution vector and other local variables
  type(N_Vector), pointer  :: sunvec_y
  type(SUNLinearSolver), pointer :: sunlinsol_LS
  type(SUNMatrix), pointer :: sunmat_A
  type(c_ptr)              :: arkode_mem           ! ARKODE memory structure
  type(c_ptr)              :: sunctx               ! SUNDIALS context structure
  integer, pointer         :: commptr              ! comm pointer to initialize sunctx

  double precision   :: t(1), dTout, tout, urms
  integer(c_long)    :: N, Ntot
  integer(c_int)     :: flag, retval, ioutput
  logical            :: outproc
  character(len=100) :: outname

  ! initialize MPI
  call MPI_Init(flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Init = ", flag
     stop
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Comm_rank = ", flag
     call MPI_Finalize(flag)
  end if

  ! Create SUNDIALS simulation context
  comm = MPI_COMM_WORLD
  commptr => comm
  retval = FSUNContext_Create(c_loc(commptr), sunctx)
  if (retval /= 0) then
    print *, "Error: FSUNContext_Create returned ",retval
    call MPI_Finalize(flag)
  end if

!   ! Configure SUNDIALS logging via the runtime API.
!   ! This requires that SUNDIALS was configured with the CMake options
!   !   SUNDIALS_LOGGING_LEVEL=n     (see logger in HeatUserData above)
!   ! where n is one of:
!   !    1 --> log only errors,
!   !    2 --> log errors + warnings,
!   !    3 --> log errors + warnings + info output
!   !    4 --> all of the above plus debug output
!   !    5 --> all of the above and even more
!   ! SUNDIALS will only log up to the max level n, but a lesser level can
!   ! be configured at runtime by only providing output files for the
!   ! desired levels. We will enable informational logging here:
!   retval = FSUNLogger_Create(c_loc(commptr), 0, logger)
!   if (retval /= 0) then
!     print *, "Error: FSUNLogger_Create returned ",retval
!     call MPI_Abort(comm, 1, ierr)
!   end if

  ! Initialize HeatUserData module
  call InitUserData()
  nx = nx_
  ny = ny_
  kx = kx_
  ky = ky_
  dx = 1.d0/(nx-1)   ! x mesh spacing
  dy = 1.d0/(ny-1)   ! x mesh spacing

  ! Set up parallel decomposition (computes local mesh sizes)
  call SetupDecomp(flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in SetupDecomp = ", flag
     call MPI_Finalize(flag)
  end if

  ! Initialize data structures
  N = nxl*nyl
  Ntot = nx*ny
  ! Create solution vector
  sunvec_y  => FN_VNew_Parallel(comm, N, Ntot, sunctx)

  ! Initial problem output
  outproc = (myid == 0)
  if (outproc) then
     write(6,*) "  "
     write(6,*) "2D Heat PDE test problem:";
     write(6,'(A,i4)') "   nprocs = " , nprocs
     write(6,'(A,i4)') "   nx = ", nx
     write(6,'(A,i4)') "   ny = ", ny
     write(6,'(A,f5.2)') "   kx = ", kx
     write(6,'(A,f5.2)') "   ky = ", ky
     write(6,'(A,es9.2)') "   rtol = ", rtol
     write(6,'(A,es9.2)') "   atol = ", atol
     write(6,'(A,i4)') "   nxl (proc 0) = ", nxl
     write(6,'(A,i4)') "   nyl (proc 0) = ", nyl
     write(6,*) "  "
  end if

  ! Set initial conditions
  call FN_VConst(0.d0, sunvec_y)

  ! initialize PCG linear solver module
  sunlinsol_LS => FSUNLinSol_PCG(sunvec_y, PCGpretype, PCGmaxl, sunctx)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in FSunLinSol_PCG = ", flag
     call MPI_Finalize(flag)
  end if

  arkode_mem = FARKStepCreate(c_null_ptr, c_funloc(farkifun), T0, sunvec_y, sunctx)
  if (.not. c_associated(arkode_mem)) then
     print *, "Error: FARKStepCreate returned NULL"
     call MPI_Finalize(flag)
  end if

  ! Initialize the scalar tolerances for the solver
  retval = FARKStepSStolerances(arkode_mem, rtol, atol)
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSStolerances = ", retval, "; halting"
     call MPI_Finalize(flag)
  end if

  ! fill in the heat source array
  do j=1,nyl
     do i=1,nxl
        h(i,j) = sin(pi*(is+i-2)*dx) * sin(2.d0*pi*(js+j-2)*dy)
     end do
  end do

  ! set integrator options
  retval = FARKStepSetNonlinConvCoef(arkode_mem, nlscoef)
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSetNonlinConvCoef = ", retval
     call MPI_Finalize(flag)
  end if
  retval = FARKStepSetPredictorMethod(arkode_mem, 1)
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSetPredictorMethod = ", retval
     call MPI_Finalize(flag)
  end if

  ! attach linear solver module to ARKLs interface

  sunmat_A => null()

  retval = FARKStepSetLinearSolver(arkode_mem, sunlinsol_LS, sunmat_A)
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSetLinearSolver = ", retval
     call MPI_Finalize(flag)
  end if

  ! Signal user-supplied preconditioner
  retval = FARKStepSetPreconditioner(arkode_mem, c_funloc(farkpset), c_funloc(farkpsol))
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSetPreconditioner = ", retval
     call MPI_Finalize(flag)
  end if

  ! specify that the problem is linearly implicit, but that Jacobian does not depend on time
  retval = FARKStepSetLinear(arkode_mem, 0)
  if (retval /= 0) then
     write(0,*) "Error in FARKStepSetLinear = ", retval
     call MPI_Finalize(flag)
  end if

  ! Each processor outputs subdomain information
  write(outname,'(16Hheat2d_subdomain,f4.3,4H.txt)') myid/1000.0
  open(100, file=outname)
  write(100,'(6(i9,1x))') nx, ny, is, ie, js, je
  close(100)

  ! Open output streams for results, access data array
  write(outname,'(6Hheat2d,f4.3,4H.txt)') myid/1000.0
  open(101, file=outname)

  call InitOutput_Disk(sunvec_y, 101, nxl, nyl)

  ! Main time-stepping loop: calls ARKode to perform the integration, then
  ! prints results.  Stops when the final time has been reached
  t = T0
  dTout = (Tf-T0)/Nt
  tout = T0+dTout
  call PRMS(sunvec_y, urms, flag)
  if (outproc) then
    write(6,*) "        t      ||u||_rms"
    write(6,*) "   ----------------------"
    write(6,'(2(2x,f10.6))') t, urms
  end if
  do ioutput=1,Nt

     retval = FARKStepEvolve(arkode_mem, tout, sunvec_y, t, ARK_NORMAL)         ! call integrator
     if (retval /= 0) then
        write(0,*) "Error in FARKStepEvolve = ", retval
        call MPI_Finalize(flag)
     end if

     call PRMS(sunvec_y, urms, flag)
     if (outproc) &
          write(6,'(2(2x,f10.6))') t, urms     ! print solution stats
     if (flag >= 0) then                       ! successful solve: update output time
        tout = min(tout + dTout, Tf)
     else                                      ! unsuccessful solve: break
        if (outproc) &
             write(0,*) "Solver failure, stopping integration"
        exit
     end if

     call OutputResults_Disk(sunvec_y, 101, nxl, nyl)

  end do
  if (outproc) then
     write(6,*) "   ----------------------"
  end if
  close(101)

  ! Print some final statistics
  if (outproc) then
     call GetFinalStats(arkode_mem, flag, retval)
     write(6,*) "  "
     write(6,*) "Final Solver Statistics:"
     write(6,'(2(A,i6),A)') "   Internal solver steps = ", nst(1), &
          " (attempted = ", nst_att(1), ")"
     write(6,'(A,i6,A,i6)') "   Total RHS evals:  Fe = ", nfe(1), ",  Fi = ", nfi(1)
     write(6,'(A,i6)') "   Total linear solver setups = ", nls(1)
     write(6,'(A,i6)') "   Total linear iterations = ", nli(1)
     write(6,'(A,i6)') "   Total number of Jacobian-vector products = ", njve(1)
     write(6,'(A,i6)') "   Total number of Preconditioner setups = ", npe(1)
     write(6,'(A,i6)') "   Total number of Preconditioner solves = ", nps(1)
     write(6,'(A,i6)') "   Total number of linear solver convergence failures = ", &
          nlcf(1)
     write(6,'(A,i6)') "   Total number of Newton iterations = ", niters(1)
     write(6,'(A,i6)') "   Total number of nonlinear solver convergence failures = ", &
          nlscf(1)
     write(6,'(A,i6)') "   Total number of error test failures = ", netf(1)
 end if

 ! Clean up and return with successful completion
 call FreeUserData(flag)               ! free user data
 call FARKStepFree(arkode_mem)
 call FN_VDestroy(sunvec_y)
 retval = FSUNLinSolFree(sunlinsol_LS) ! free linear solver
 if (retval /= 0) then
    print *, "Error: FSUNLinSolFree returned ", retval
    call MPI_Finalize(flag)
 end if
 call MPI_Barrier(comm, flag)
 call MPI_Finalize(flag)               ! Finalize MPI

end program driver
!-----------------------------------------------------------------
