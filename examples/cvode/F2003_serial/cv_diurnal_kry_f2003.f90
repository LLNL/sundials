! ------------------------------------------------------------------
! Programmer(s): Daniel M. Margolis @ SMU
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
! FCVODE Example Problem: 2D kinetics-transport,
! precond. Krylov solver.
!
! An ODE system is generated from the following 2-species diurnal
! kinetics advection-diffusion PDE system in 2 space dimensions:
!
! dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
!                       + Ri(c1,c2,t)      for i = 1,2,   where
! R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
! R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
! Kv(y) = Kv0*exp(y/5) ,
! Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
! vary diurnally.
!
! The problem is posed on the square
! 0 <= x <= 20,    30 <= y <= 50   (all in km),
! with homogeneous Neumann boundary conditions, and for time t in
! 0 <= t <= 86400 sec (1 day).
!
! The PDE system is treated by central differences on a uniform
! 10 x 10 mesh, with simple polynomial initial profiles.
! The problem is solved with CVODE, with the BDF/GMRES method and
! and the block-diagonal part of the Jacobian as a left preconditioner.
!
! Note that this problem should only work with SUNDIALS configured
! to use 'realtype' as 'double' and 'sunindextype' as '64bit'
!
! The second and third dimensions of U here must match the values of
! MX and MY, for consistency with the output statements below.
! ------------------------------------------------------------------

module diurnal_mod

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! setup and number of equations
    integer(c_int), parameter  :: mx = 10, my = 10
    integer(c_long), parameter :: mm = mx*my
    integer(c_long), parameter :: neq = 2*mm

    ! ODE constant parameters
    double precision, parameter :: Kh = 4.0d-6
    double precision, parameter :: Vel = 0.001d0
    double precision, parameter :: Kv0 = 1.0d-8
    double precision, parameter :: q1 = 1.63d-16
    double precision, parameter :: q2 = 4.66d-16
    double precision, parameter :: c3 = 3.7d16
    double precision, parameter :: pi = 3.1415926535898d0
    double precision, parameter :: halft = 4.32d4
    double precision, parameter :: om = pi/halft
    double precision, parameter :: dx = 20.0d0/(mx - 1.0d0)
    double precision, parameter :: dy = 20.0d0/(my - 1.0d0)
    double precision, parameter :: hdco = Kh/(dx**2)
    double precision, parameter :: haco = Vel/(2.0d0*dx)
    double precision, parameter :: vdco = (1.0d0/(dy**2))*Kv0
    double precision, parameter :: a3 = 22.62d0
    double precision, parameter :: a4 = 7.601d0

    ! Solving assistance fixed parameters
    double precision, parameter :: twohr = 7200.0D0
    double precision, parameter :: rtol = 1.0d-5
    double precision, parameter :: floor = 100.0d0
    double precision, parameter :: delt = 0.0d0
    double precision, parameter :: atol = rtol*floor
    integer(c_int), parameter   :: Jpretype = 1
    integer(c_int), parameter   :: iGStype = 1
    integer(c_int), parameter   :: maxL = 0
    integer(c_long), parameter  :: mxsteps = 1000
    double precision :: p_p(2,2,mm)

    ! ODE non-constant parameters
    double precision :: q3
    double precision :: q4
    double precision :: c1
    double precision :: c2
    integer(c_long)  :: jx, jy

    contains

    ! ----------------------------------------------------------------
    ! RhsFn provides the right hand side implicit function for the
    ! ODE: dy1/dt = f1(t,y1,y2,y3)
    !      dy2/dt = f2(t,y1,y2,y3)
    !      dy3/dt = f3(t,y1,y2,y3)
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer(c_int) function RhsFn(tn, sunvec_u, sunvec_f, user_data) &
         result(ierr) bind(C,name='RhsFn')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      double precision, value :: tn    ! current time
      type(N_Vector) :: sunvec_u       ! solution N_Vector
      type(N_Vector) :: sunvec_f       ! rhs N_Vector
      type(c_ptr) :: user_data         ! user-defined data

      ! local data
      integer(c_int)   :: ileft, iright
      integer(c_long)  :: idx0, idn, iup, idx
      double precision :: c1dn, c2dn, c1up, c2up, c1lt, c2lt
      double precision :: c1rt, c2rt, cydn, cyup, hord1, hord2, horad1
      double precision :: horad2, qq1, qq2, qq3, qq4, rkin1, rkin2, s
      double precision :: vertd1, vertd2, ydn, yup

      ! pointers to data in SUNDIALS vectors
      double precision, pointer, dimension(2,mm) :: uvec(:,:)
      double precision, pointer, dimension(2,mm) :: fvec(:,:)

      !======= Internals ============

      ! get data arrays from SUNDIALS vectors
      uvec(1:2,1:mm) => FN_VGetArrayPointer(sunvec_u)
      fvec(1:2,1:mm) => FN_VGetArrayPointer(sunvec_f)

      ! Set diurnal rate coefficients.
      s = sin(om * tn)
      if (s > 0.0d0) then
         q3 = exp(-a3 / s)
         q4 = exp(-a4 / s)
      else
         q3 = 0.0d0
         q4 = 0.0d0
      end if

      ! Loop over all grid points.
      do jy = 1, my
         ydn = 30.0d0 + (jy - 1.5d0) * dy
         yup = ydn + dy
         cydn = vdco * exp(0.2d0 * ydn)
         cyup = vdco * exp(0.2d0 * yup)
         idx0 = (jy - 1) * mx
         idn = -mx
         if (jy == 1) idn = mx
         iup = mx
         if (jy == my) iup = -mx
         do jx = 1, mx
            idx = idx0 + jx
            c1 = uvec(1,idx)
            c2 = uvec(2,idx)
            ! Set kinetic rate terms.
            qq1 = q1 * c1 * c3
            qq2 = q2 * c1 * c2
            qq3 = q3 * c3
            qq4 = q4 * c2
            rkin1 = -qq1 - qq2 + 2.0d0 * qq3 + qq4
            rkin2 = qq1 - qq2 - qq4
            ! Set vertical diffusion terms.
            c1dn = uvec(1,idx + idn)
            c2dn = uvec(2,idx + idn)
            c1up = uvec(1,idx + iup)
            c2up = uvec(2,idx + iup)
            vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn)
            vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn)
            ! Set horizontal diffusion and advection terms.
            ileft = -1
            if (jx == 1) ileft = 1
            iright = 1
            if (jx == mx) iright = -1
            c1lt = uvec(1,idx + ileft)
            c2lt = uvec(2,idx + ileft)
            c1rt = uvec(1,idx + iright)
            c2rt = uvec(2,idx + iright)
            hord1 = hdco * (c1rt - 2.0d0 * c1 + c1lt)
            hord2 = hdco * (c2rt - 2.0d0 * c2 + c2lt)
            horad1 = haco * (c1rt - c1lt)
            horad2 = haco * (c2rt - c2lt)
            ! load all terms into fvec.
            fvec(1,idx) = vertd1 + hord1 + horad1 + rkin1
            fvec(2,idx) = vertd2 + hord2 + horad2 + rkin2
         end do
      end do

      ! return success
      ierr = 0
      return

    end function RhsFn
    
    integer(c_int) function PreSet(t, sunvec_u, sunvec_f, jok, jcur, gamma, user_data) &
            result(ierr) bind(C,name='PreSet')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      double precision,  value :: t
      type(N_Vector)   :: sunvec_u       ! solution N_Vector
      type(N_Vector)   :: sunvec_f       ! f-value N_Vector
      integer(c_int), value :: jok
      integer(c_int)   :: jcur
      double precision, value :: gamma
      type(c_ptr) :: user_data

      ! temporary variables
      double precision, pointer, dimension(2,mm)   :: u(:,:)
      double precision :: p_bd(2,2,mm)
      u(1:2,1:mm) => FN_VGetArrayPointer(sunvec_u)

      ! if needed, recompute bd
      if (jok == 1) then
      ! jok = 1. reuse saved bd
         jcur = 0
      else    
      ! jok = 0. compute diagonal jacobian blocks.
      ! (using q4 value computed on last fcvfun call).
         call Prec_Jac(mx, my, u, p_bd, q1, q2, q3, q4, &
                       c3, dy, hdco, vdco, ierr)
         jcur = 1
      endif

      ! copy bd to p and scale by -gamma
      p_p = -gamma * p_bd

      ! Perform LU decomposition
      call Prec_LU(mm, p_p, ierr)

      ! return success
      ierr = 0
      return

    end function PreSet

    integer(c_int) function PreSolve(t, sunvec_u, sunvec_f, sunvec_r, sunvec_z, &
            gamma, delta, lr, user_data) result(ierr) bind(C,name='PreSolve')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      double precision, value :: t
      type(N_Vector) :: sunvec_u       ! solution N_Vector
      type(N_Vector) :: sunvec_f       ! f-value N_Vector
      type(N_Vector) :: sunvec_r       ! rhs N_Vector
      type(N_Vector) :: sunvec_z       ! output N_Vector
      double precision :: gamma, delta
      integer(c_int) :: lr
      type(c_ptr)  :: user_data

      ! temporary variables
      double precision, pointer, dimension(2,mm) :: z(:,:), r(:,:)
      z(1:2,1:mm) => FN_VGetArrayPointer(sunvec_z)
      r(1:2,1:mm) => FN_VGetArrayPointer(sunvec_r)

      ! copy rhs into z
      z = r

      ! solve the block-diagonal system px = r using lu factors stored in p
      ! and pivot data in ipp, and return the solution in z.
      call Prec_Sol(mm, p_p, z)

      ! return success
      ierr = 0
      return
      
    end function PreSolve

    subroutine Prec_Jac(mmx, mmy, u, bd, qq1, qq2, qq3, qq4, cc3, &
                        ddy, hhdco, vvdco, ierr)
      
      implicit none

      integer(c_int), intent(in)    :: mmx, mmy
      double precision, intent(in)  :: u(2,*)
      double precision, intent(out) :: bd(2,2,*)
      double precision, intent(in)  :: qq1, qq2, qq3, qq4, cc3, ddy, hhdco, vvdco
      integer(c_int), intent(out)   :: ierr

      ! local variables
      integer(c_int)   :: idx, idx0
      double precision :: cydn, cyup, diag, ydn, yup

      do jy = 1, mmy
         ydn = 30.0d0 + (jy - 1.5d0) * ddy
         yup = ydn + ddy
         cydn = vvdco * exp(0.2d0 * ydn)
         cyup = vvdco * exp(0.2d0 * yup)
         diag = -(cydn + cyup + 2.0d0 * hhdco)
         idx0 = (jy - 1) * mmx
         do jx = 1, mmx
            idx = idx0 + jx
            c1 = u(1,idx)
            c2 = u(2,idx)
            bd(1,1,idx) = (-qq1 * cc3 - qq2 * c2) + diag
            bd(1,2,idx) = -qq2 * c1 + qq4
            bd(2,1,idx) =  qq1 * cc3 - qq2 * c2
            bd(2,2,idx) = (-qq2 * c1 - qq4) + diag
         end do
      end do
      
      ! return success
      ierr = 0
      return
   
    end subroutine Prec_Jac
    
    subroutine Prec_LU(mmm, p, ierr)

      implicit none

      integer(c_int), intent(out) :: ierr
      integer(c_long), intent(in) :: mmm
      double precision, intent(inout) :: p(2,2,mmm)

      ! local variable
      integer(c_long)  :: i
      double precision :: p11, p12, p21, p22, det

      ! add identity matrix and do lu decompositions on blocks, in place.
      do i = 1,mmm
         p11 = p(1,1,i) + 1.0d0
         p22 = p(2,2,i) + 1.0d0
         p12 = p(1,2,i)
         p21 = p(1,2,i)
         det = p11*p22 - p12*p21
         if (det == 0.d0) return

         p(1,1,i) = p22/det
         p(2,2,i) = p11/det
         p(1,2,i) = -p21/det
         p(2,1,i) = -p12/det
      end do

      return
   
    end subroutine Prec_LU

    subroutine Prec_Sol(mmm, p, z)
    
      implicit none

      integer(c_long), intent(in)   :: mmm
      double precision, dimension(2,2,mmm), intent(inout) :: p(:,:,:)
      double precision, dimension(2,mmm), intent(inout) :: z(:,:)
      
      ! local variable
      integer(c_long)  :: i
      double precision :: z1, z2


      do i = 1,mm
         z1 = z(1,i)
         z2 = z(2,i)
         z(1,i) = p(1,1,i) * z1 + p(1,2,i) * z2
         z(2,i) = p(2,1,i) * z1 + p(2,2,i) * z2
      end do

      return

    end subroutine Prec_Sol
    
  end module diurnal_mod

  program main

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    use fcvode_mod                 ! Fortran interface to the ARKode module
    use fsundials_futils_mod       ! Fortran utilities
    use fsundials_nvector_mod      ! Fortran interface to the generic N_Vector
    use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
    use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
    use fnvector_serial_mod        ! Fortran interface to serial N_Vector
    use fsunlinsol_spgmr_mod       ! Fortran interface to spgmr SUNLinearSolver
    use fsundials_context_mod      ! Fortran interface to SUNContext
    use diurnal_mod                ! ODE functions

    !======= Declarations =========
    implicit none

    ! local variables
    type(c_ptr)      :: ctx        ! SUNDIALS context for the simulation
    double precision :: tstart     ! initial time
    double precision :: tout       ! output time
    double precision :: tcur(1)    ! current time
    double precision :: cx, cy     ! initialization variables
    integer(c_int)   :: ierr       ! error flag from C functions
    integer(c_long)  :: outstep    ! output step
    double precision :: x, y       ! initialization index variables

    type(N_Vector),        pointer :: sunvec_u      ! sundials vector
    type(N_Vector),        pointer :: sunvec_f      ! sundials vector
    type(SUNLinearSolver), pointer :: sunls         ! sundials linear solver
    type(SUNMatrix),       pointer :: sunmat_A      ! sundials matrix (empty)
    type(c_ptr)                    :: cvode_mem     ! CVODE memory
    double precision, pointer, dimension(2,mx,my) :: uvec(:,:,:) ! underlying vector

    ! output statistic variables
    integer(c_long)  :: lnst(1)
    double precision :: lh(1)

    !======= Internals ============

    ! create the SUNDIALS context
    ierr = FSUNContext_Create(c_null_ptr, ctx)

    ! initialize ODE
    tstart = 0.0d0
    tcur   = tstart

    ! create SUNDIALS N_Vector
    sunvec_u => FN_VNew_Serial(neq, ctx)
    if (.not. associated(sunvec_u)) then
       print *, 'ERROR: sunvec = NULL'
       stop 1
    end if

    sunvec_f => FN_VNew_Serial(neq, ctx)
    if (.not. associated(sunvec_f)) then
       print *, 'ERROR: sunvec = NULL'
       stop 1
    end if

    uvec(1:2,1:mx,1:my) => FN_VGetArrayPointer(sunvec_u)

    ! initialize and fill initial condition vector
    do jy = 1,my
       y  = 30.0d0 + (jy - 1.0d0) * dy
       cy = (0.1d0 * (y - 40.0d0))**2
       cy = 1.0d0 - cy + 0.5d0 * cy**2
       do jx = 1,mx
          x = (jx - 1.0d0) * dx
          cx = (0.1d0 * (x - 10.0d0))**2
          cx = 1.0d0 - cx + 0.5d0 * cx**2
          uvec(1,jx,jy) = 1.0d6 * cx * cy
          uvec(2,jx,jy) = 1.0d12 * cx * cy
       end do
    end do

    ! create and initialize CVode memory
    cvode_mem = FCVodeCreate(CV_BDF, ctx)
    if (.not. c_associated(cvode_mem)) print *,'ERROR: cvode_mem = NULL'

    ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_u)
    if (ierr /= 0) then
       print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeSetMaxNumSteps(cvode_mem, mxsteps)
    if (ierr /= 0) then
       print *, 'Error in FCVodeSetMaxNumSteps, ierr = ', ierr, '; halting'
       stop 1
    end if
    
    ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
    if (ierr /= 0) then
       print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
       stop 1
    end if

    ! Tell CVODE to use a SPGMR linear solver.
    sunls => FSUNLinSol_SPGMR(sunvec_u, Jpretype, maxL, ctx)
    if (.not. associated(sunls)) then
       print *, 'ERROR: sunls = NULL'
       stop 1
    end if

    ierr = FSUNLinSol_SPGMRSetGSType(sunls, iGStype)
    if (ierr /= 0) then
       print *, 'Error in FCVodeSetLinearSolver'
       stop 1
    end if

    ! Attach the linear solver (with NULL SUNMatrix object)
    sunmat_A => null()
    ierr = FCVodeSetLinearSolver(cvode_mem, sunls, sunmat_A)
    if (ierr /= 0) then
       print *, 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeSetPreconditioner(cvode_mem, c_funloc(PreSet), c_funloc(PreSolve))
    if (ierr /= 0) then
      print *, 'Error in FCVodeSetPreconditioner, ierr = ', ierr, '; halting'
      stop 1
   end if

    ! Start time stepping
    print *, '   '
    print *, 'Finished initialization, starting time steps'
    print *, '   '
    print *, '      t          c1  (bottom left      middle       top right)  | lnst  lh'
    print *, '      t          c2  (bottom left      middle       top right)  | lnst  lh'
    print *, ' -----------------------------------------------------------------------------------'
    tout = twohr
    do outstep = 1,12

       ! call CVode
       ierr = FCVode(cvode_mem, tout, sunvec_u, tcur, CV_NORMAL)
       if (ierr /= 0) then
          print *, 'Error in FCVode, ierr = ', ierr, '; halting'
          stop 1
       end if

       ierr = FCVodeGetNumSteps(cvode_mem, lnst)
       if (ierr /= 0) then
          print *, 'Error in FCVodeGetNumSteps, ierr = ', ierr, '; halting'
          stop 1
       end if

       ierr = FCVodeGetCurrentStep(cvode_mem, lh)
       if (ierr /= 0) then
          print *, 'Error in FCVodeGetCurrentStep, ierr = ', ierr, '; halting'
          stop 1
       end if

       ! print current solution and output statistics
       print '(2x,4(es14.6,2x),i5,es14.6)', tcur, uvec(1,1,1), uvec(1,5,5), uvec(1,10,10), lnst, lh
       print '(18x,3(es14.6,2x))', uvec(2,1,1), uvec(2,5,5), uvec(2,10,10)

       ! update tout
       tout = tout + twohr

    end do
    print *, ' -----------------------------------------------------------------------------------'

    ! diagnostics output
    call CVodeStats(cvode_mem)

    ! clean up
    call FCVodeFree(cvode_mem)
    call FN_VDestroy(sunvec_u)
    call FN_VDestroy(sunvec_f)
    ierr = FSUNLinSolFree(sunls)
    ierr = FSUNContext_Free(ctx)

  end program main

  ! ----------------------------------------------------------------
  ! CVodeStats
  !
  ! Print CVODE statstics to stdandard out
  ! ----------------------------------------------------------------
  subroutine CVodeStats(cvode_mem)

    !======= Inclusions ===========
    use iso_c_binding
    use fcvode_mod

    !======= Declarations =========
    implicit none

    type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

    integer(c_int)  :: ierr          ! error flag

    integer(c_long) :: nsteps(1)     ! num steps
    integer(c_long) :: nfe(1)        ! num function evals
   !  integer(c_long) :: npsetups(1)   ! num preconditioner setups
    integer(c_long) :: netfails(1)   ! num error test fails
    integer(c_long) :: npe(1)        ! num preconditioner evals
    integer(c_long) :: nps(1)        ! num preconditioner solves
    integer(c_long) :: nniters(1)    ! nonlinear solver iterations
    integer(c_long) :: nliters(1)    ! linear solver iterations
   !  integer(c_long) :: nmcf(1)       ! num mass convergence failures
    integer(c_long) :: ncf(1)        ! num convergence failures nonlinear
    integer(c_long) :: ncfl(1)       ! num convergence failures linear
    integer(c_long) :: nncfails(1)   ! nonlinear solver fails
    integer(c_long) :: lenrw(1)      ! main solver real/int workspace size
    integer(c_long) :: leniw(1)
    integer(c_long) :: lenrwls(1)    ! linear solver real/int workspace size
    integer(c_long) :: leniwls(1)
    integer(c_long) :: nfebp(1)      ! num f evaluations
    double precision :: avdim(1)     ! avg Krylov subspace dim (NLI/NNI)
    integer(c_long) :: lenrwbp(1)
    integer(c_long) :: leniwbp(1)

    !======= Internals ============

    ierr = FCVodeGetNumSteps(cvode_mem, nsteps)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumSteps, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumRhsEvals(cvode_mem, nfe)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumRhsEvals, ierr = ', ierr, '; halting'
       stop 1
    end if

   !  ierr = FCVodeGetNumMassSetups(cvode_mem, npsetups)
   !  if (ierr /= 0) then
   !     print *, 'Error in FCVodeGetNumMassSetups, ierr = ', ierr, '; halting'
   !     stop 1
   !  end if

    ierr = FCVodeGetNumErrTestFails(cvode_mem, netfails)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumErrTestFails, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumPrecEvals(cvode_mem, npe)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumPrecEvals, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumPrecSolves(cvode_mem, nps)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumPrecSolves, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumNonlinSolvIters(cvode_mem, nniters)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumNonlinSolvIters, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumLinIters(cvode_mem, nliters)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumLinIters, ierr = ', ierr, '; halting'
       stop 1
    end if

    avdim = dble(nliters)/dble(nniters)

   !  ierr = FCVodeGetNumMassConvFails(cvode_mem, nmcf)
   !  if (ierr /= 0) then
   !     print *, 'Error in FCVodeGetNumMassConvFails, ierr = ', ierr, '; halting'
   !     stop 1
   !  end if

    ierr = FCVodeGetNumLinConvFails(cvode_mem, ncfl)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumLinConvFails, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetNumNonlinSolvConvFails(cvode_mem, ncf)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetNumNonlinSolvConvFails, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetWorkSpace(cvode_mem, lenrw, leniw)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetWorkSpace, ierr = ', ierr, '; halting'
       stop 1
    end if

    ierr = FCVodeGetLinWorkSpace(cvode_mem, lenrwls, leniwls)
    if (ierr /= 0) then
       print *, 'Error in FCVodeGetLinWorkSpace, ierr = ', ierr, '; halting'
       stop 1
    end if

    print *, ' '
    print *, ' General Solver Stats:'
    print '(4x,A,i9)'       ,'Total internal steps taken      =',nsteps
    print '(4x,A,i9)'       ,'Total rhs function call         =',nfe
   !  print '(4x,A,i9)'       ,'Num lin solver setup calls    =',npsetups
    print '(4x,A,i9)'       ,'Total num preconditioner evals  =',npe
    print '(4x,A,i9)'       ,'Total num preconditioner solves =',nps
    print '(4x,A,i9)'       ,'Num error test failures         =',netfails
    print '(4x,A,i9)'       ,'Num nonlinear solver iters      =',nniters
    print '(4x,A,i9)'       ,'Num linear solver iters         =',nliters
    print '(4x,A,es14.6)'   ,'Avg Krylov subspace dim         =',avdim
    print '(4x,A,i9)'       ,'Num nonlinear solver fails      =',ncf
    print '(4x,A,i9)'       ,'Num linear solver fails         =',ncfl
   !  print '(4x,A,i9)'       ,'Num mass solver fails         =',nmcf
    print '(4x,A,2(i9,3x))' ,'main solver real/int workspace sizes   =',lenrw,leniw
    print '(4x,A,2(i9,3x))' ,'linear solver real/int workspace sizes =',lenrwls,leniwls
    print *, ' '

    return

  end subroutine CVodeStats
