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
! Example problem:
!
! The following test simulates a brusselator problem from chemical
! kinetics.  This is a PDE system with 3 components, Y = [u,v,w],
! satisfying the equations,
!    du/dt = du*u_xx + a - (w+1)*u + v*u^2
!    dv/dt = dv*v_xx + w*u - v*u^2
!    dw/dt = dw*w_xx + (b-w)/ep - w*u
! for t in the interval [0, 10], x in [0, 10], with initial
! conditions
!    u(0,x) =  a  + 0.1*sin(pi*x),
!    v(0,x) = b/a + 0.1*sin(pi*x),
!    w(0,x) =  b  + 0.1*sin(pi*x),
! and with stationary boundary conditions, i.e.
!    u_t(t,0) = u_t(t,1) = 0,
!    v_t(t,0) = v_t(t,1) = 0,
!    w_t(t,0) = w_t(t,1) = 0.
!
! Here, we use a piecewise linear Galerkin finite element
! discretization in space, where all element-wise integrals are
! computed using 3-node Gaussian quadrature (since we will have
! quartic polynomials in the reaction terms for the u_t and v_t
! equations (including the test function)).  The time derivative
! terms in this system will include a mass matrix, giving rise to
! an ODE system of the form
!      M y_t = L y + R(y),
! where M is the 3x3 block mass matrix for each component, L is
! the 3x3 block Laplace operator for each component, and R(y) is
! comprised of the nonlinear reaction terms for each component.
! Since it it highly inefficient to rewrite this system as
!      y_t = M^{-1}(L y + R(y)),
! we solve this system using FARKODE, with a user-supplied mass
! matrix.  We therefore provide functions to evaluate the ODE RHS
!    f(t,y) = L y + R(y),
! its Jacobian
!    J(t,y) = L + dR/dy,
! and the mass matrix, M.
!
! We use N=201 spatial nodes, with parameters
!    a=0.6,  b=2.0,  du=0.025,  dv=0.025,  dw=0.025,  ep=1.d-5
!
! This program solves the problem with the DIRK method, using a
! Newton iteration with the SUNKLU sparse linear solvers for both
! the system and mass matrices.  These matrices are stored in
! compressed-sparse-row format.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
! ------------------------------------------------------------------

module UserData

   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
 
   !======= Declarations =========
   implicit none
 
   ! number of equations
   integer(c_long), parameter :: neqreal = 3
 
   ! ODE parameters
   integer(c_long), parameter  :: N = 201           ! number of intervals
   integer(c_long), parameter  :: neq = neqreal*N   ! set overall problem size
   double precision, parameter :: a = 0.6d0         ! constant forcing on u
   double precision, parameter :: b = 2.d0          ! steady-state value of w
   double precision, parameter :: du = 2.5d-2       ! diffusion coeff for u
   double precision, parameter :: dv = 2.5d-2       ! diffusion coeff for v
   double precision, parameter :: dw = 2.5d-2       ! diffusion coeff for w
   double precision, parameter :: ep = 1.d-5        ! stiffness parameter
   double precision, dimension(N) :: x              ! mesh node locations
 
 contains
 
   ! function that maps 2D data into 1D address space
   ! (0-based since CSR matrix will be sent to C solver)
   integer(c_long) function idx(ix, ivar)
     integer(c_long) :: ivar, ix
     idx = neqreal*(ix - 1) + ivar - 1
   end function idx

end module UserData
   
   
! finite element basis functions
module FEM

contains

   ! left/right basis functions
   double precision function ChiL(xl,xr,x)
      double precision :: xl, xr, x
      ChiL = (xr-x)/(xr-xl)
   end function ChiL

   double precision function ChiR(xl,xr,x)
      double precision :: xl, xr, x
      ChiR = (x-xl)/(xr-xl)
   end function ChiR

   ! derivatives of left/right basis functions
   double precision function ChiL_x(xl,xr)
      double precision :: xl, xr
      ChiL_x = 1.d0/(xl-xr)
   end function ChiL_X

   double precision function ChiR_x(xl,xr)
      double precision :: xl, xr
      ChiR_x = 1.d0/(xr-xl)
   end function ChiR_x

   ! FEM output evaluation routines: value and derivative
   double precision function Eval(ul,ur,xl,xr,x)
      double precision :: ul, ur, xl, xr, x
      Eval = ul*ChiL(xl,xr,x) + ur*ChiR(xl,xr,x)
   end function Eval

   double precision function Eval_x(ul,ur,xl,xr)
      double precision :: ul, ur, xl, xr
      Eval_x = ul*ChiL_x(xl,xr) + ur*ChiR_x(xl,xr)
   end function Eval_x

end module FEM

! quadrature data
module Quadrature

contains

  ! nodes
  double precision function X1(xl,xr)
     double precision :: xl, xr
     X1 = 0.5d0*(xl+xr) - 0.5d0*(xr-xl)*0.774596669241483377035853079956d0
  end function X1

  double precision function X2(xl,xr)
     double precision :: xl, xr
     X2 = 0.5d0*(xl+xr)
  end function X2

  double precision function X3(xl,xr)
     double precision :: xl, xr
     X3 = 0.5d0*(xl+xr) + 0.5d0*(xr-xl)*0.774596669241483377035853079956d0
  end function X3

  ! quadrature
  double precision function Quad(f1,f2,f3,xl,xr)
    real(c_double) :: f1, f2, f3, xl, xr
    real(c_double), parameter :: wt1=0.55555555555555555555555555555556d0
    real(c_double), parameter :: wt2=0.88888888888888888888888888888889d0
    real(c_double), parameter :: wt3=0.55555555555555555555555555555556d0
    Quad = 0.5d0*(xr-xl)*(wt1*f1 + wt2*f2 + wt3*f3)
  end function Quad

end module Quadrature

module ode_mod

   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
   use UserData

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
   use FEM
   use Quadrature
   use fsundials_nvector_mod

   !======= Declarations =========
   implicit none

   ! calling variables
   real(c_double), value :: tn                  ! current time
   type(N_Vector) :: sunvec_y                   ! solution N_Vector
   type(N_Vector) :: sunvec_f                   ! rhs N_Vector
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
   use FEM
   use Quadrature   
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
!-----------------------------------------------------------------




!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
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
   use UserData                   ! Declarations and indexing
 
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
   integer(c_int) :: ordering, sparsetype     ! AMD and CSR types respectively
   integer(c_long):: mxsteps                  ! max num steps
   real(c_double) :: pi, h, z                 ! constants and variables to help with mesh

   ! real(c_double), parameter :: nlscoef = 1.d-2  ! non-linear solver coefficient
   ! integer(c_int), parameter :: order = 3        ! method order
 
   type(N_Vector),        pointer :: sunvec_y    ! sundials vector
   type(N_Vector),        pointer :: sunvec_u    ! sundials vector
   type(N_Vector),        pointer :: sunvec_v    ! sundials vector
   type(N_Vector),        pointer :: sunvec_w    ! sundials vector
   type(SUNMatrix),       pointer :: sunmat_A    ! sundials (linsol) matrix
   type(SUNMatrix),       pointer :: sunmat_M    ! sundials (mass) matrix
   type(SUNLinearSolver), pointer :: sunls       ! sundials linear solver
   type(c_ptr)                    :: arkode_mem  ! ARKODE memory
   real(c_double),        pointer :: yvec(:,:)   ! underlying vector y
   real(c_double),        pointer :: umask(:,:)  ! underlying vector u
   real(c_double),        pointer :: vmask(:,:)  ! underlying vector v
   real(c_double),        pointer :: wmask(:,:)  ! underlying vector w
 
   !======= Internals ============
 
   ! create the SUNDIALS context
   ierr = FSUNContext_Create(c_null_ptr, ctx)
 
   ! initialize ODE
   tstart = 0.0d0
   tend   = 10.0d0
   tcur   = tstart
   tout   = tstart
   dtout  = (tend - tstart)/10.d0
   nout   = ceiling(tend/dtout)
 
   ! create and assign SUNDIALS N_Vector
   sunvec_y => FN_VNew_Serial(neq, ctx)
   if (.not. associated(sunvec_y)) then
      print *, 'ERROR: sunvec = NULL'
      stop 1
   end if

   sunvec_u => FN_VNew_Serial(neq, ctx)
   if (.not. associated(sunvec_u)) then
      print *, 'ERROR: sunvec = NULL'
      stop 1
   end if

   sunvec_v => FN_VNew_Serial(neq, ctx)
   if (.not. associated(sunvec_v)) then
      print *, 'ERROR: sunvec = NULL'
      stop 1
   end if

   sunvec_w => FN_VNew_Serial(neq, ctx)
   if (.not. associated(sunvec_w)) then
      print *, 'ERROR: sunvec = NULL'
      stop 1
   end if

   yvec  => FN_VGetArrayPointer(sunvec_y)
   umask => FN_VGetArrayPointer(sunvec_u)
   vmask => FN_VGetArrayPointer(sunvec_v)
   wmask => FN_VGetArrayPointer(sunvec_w)
 
   ! set up spatial mesh; this [arbitrarily] clusters
   ! more intervals near the end points of the interval
   pi = 4.d0*atan(1.d0)
   h = 10.d0/(N - 1)
   do i=1,N
      z = -5.d0 + h*(i - 1)
      x(i) = 0.5d0/atan(5.d0)*atan(z) + 0.5d0
   end do

   ! output mesh to disk
   open(200, file='bruss_FEM_mesh.txt')
   do i=1,N
      write(200,*) x(i)
   end do
   close(200)

   ! set initial conditions
   do i=1,N
      y(1,i) =  a  + 0.1d0*sin(pi*x(i))   ! u0
      y(2,i) = b/a + 0.1d0*sin(pi*x(i))   ! v0
      y(3,i) =  b  + 0.1d0*sin(pi*x(i))   ! w0
   end do

   ! set mask values for each solution component
   umask = 0.d0
   vmask = 0.d0
   wmask = 0.d0
   do i=1,N
      umask(1,i) = 1.d0
      vmask(2,i) = 1.d0
      wmask(3,i) = 1.d0
   end do

   ! create ARKStep memory
   arkode_mem = FARKStepCreate(c_funloc(ExpRhsFn), c_funloc(ImpRhsFn), tstart, sunvec_y, ctx)
   if (.not. c_associated(arkode_mem)) print *,'ERROR: arkode_mem = NULL'
 
   ! Tell ARKODE to use a sparse linear solver.
   nnz = 15*neq
   ordering = 0
   sparsetype = 1
   sunmat_A => FSUNSparseMatrix(neq, neq, nnz, sparsetype, ctx)
   if (.not. associated(sunmat_A)) then
      print *, 'ERROR: sunmat_A = NULL'
      stop 1
   end if

   sunmat_M => FSUNSparseMatrix(neq, neq, nnz, sparsetype, ctx)
   if (.not. associated(sunmat_M)) then
      print *, 'ERROR: sunmat_M = NULL'
      stop 1
   end if
 
   sunls => FSUNLinSol_KLU(sunvec_y, sunmat_A, ctx)
   if (.not. associated(sunls)) then
      print *, 'ERROR: sunls = NULL'
      stop 1
   end if
 
   ierr = FARKStepSetMassLinearSolver(arkode_mem, sunls, sunmat_M)
   if (ierr /= 0) then
     print *, 'Error in FARKStepSetLinearSolver'
     stop 1
   end if
 
   ierr = FARKStepSetJacFn(arkode_mem, c_funloc(Jac))
   if (ierr /= 0) then
     print *, 'Error in FARKStepSetJacFn'
     stop 1
   end if

   ierr = FARKStepSetMassFn(arkode_mem, c_funloc(Mass))
   if (ierr /= 0) then
     print *, 'Error in FARKStepSetMassFn'
     stop 1
   end if
 
   ! set relative and absolute tolerances
   rtol = 1.0d-6
   atol = 1.0d-11
 
   ierr = FARKStepSStolerances(arkode_mem, rtol, atol)
   if (ierr /= 0) then
      print *, 'Error in FARKStepSStolerances, ierr = ', ierr, '; halting'
      stop 1
   end if
 
   ierr = FARKStepSetNonlinConvCoef(arkode_mem, nlscoef)
   if (ierr /= 0) then
     print *, 'Error in FARKStepSetNonlinConvCoef'
     stop 1
   end if
 
   mxsteps = 1000
   ierr = FARKStepSetMaxNumSteps(arkode_mem, mxsteps)
   if (ierr /= 0) then
     print *, 'Error in FARKStepSetNonlinConvCoef'
     stop 1
   end if
 
   imethod = 0
   idefault = 1
   pq = 0
   adapt_params = 0.d0
   ierr = FARKStepSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params)
   if (ierr /= 0) then
      print *, 'Error in FARKStepSetAdaptivityMethod, ierr = ', ierr, '; halting'
      stop 1
   end if
 
   ! Open output stream for results
   open(501, file='bruss_FEM_u.txt')
   open(502, file='bruss_FEM_v.txt')
   open(503, file='bruss_FEM_w.txt')
 
   ! output initial condition to disk
   write(501,*) ( y(1,i), i=1,N )
   write(502,*) ( y(2,i), i=1,N )
   write(503,*) ( y(3,i), i=1,N )

   ! output solver parameters to screen
   call FARKStepWriteParameters(arkode_mem, stdout)
 
   ! Start time stepping
   print *, '   '
   print *, 'Finished initialization, starting time steps'
   print *, '   '
   print *, '        t         ||u||_rms    ||v||_rms    ||w||_rms'
   print *, '  ----------------------------------------------------'
   print '(3x,4(es12.5,1x))', tcur, sqrt(sum(y*y*umask)/N), &
       sqrt(sum(y*y*vmask)/N), sqrt(sum(y*y*wmask)/N)   
   do outstep = 1,nout
 
      ! call ARKStep
      tout = min(tout + dtout, tend)
      ierr = FARKStepEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
      if (ierr < 0) then
         print *, 'Error in FARKStepEvolve, ierr = ', ierr, '; halting'
         stop 1
      endif
 
      ! output current solution information
      print '(3x,4(es12.5,1x))', Tcur, sqrt(sum(y*y*umask)/N), &
          sqrt(sum(y*y*vmask)/N), sqrt(sum(y*y*wmask)/N)
      write(501,*) ( y(1,i), i=1,N )
      write(502,*) ( y(2,i), i=1,N )
      write(503,*) ( y(3,i), i=1,N )

   end do
   print *, ' ----------------------------------------------------'
   close(501)
   close(502)
   close(503)
 
   ! diagnostics output
   call ARKStepStats(arkode_mem)
 
   ! clean up
   call FARKStepFree(arkode_mem)
   call FN_VDestroy(sunvec_y)
   call FN_VDestroy(sunvec_u)
   call FN_VDestroy(sunvec_v)
   call FN_VDestroy(sunvec_w)
   call FSUNMatDestroy(sunmat_A)
   call FSUNMatDestroy(sunmat_M)
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
 