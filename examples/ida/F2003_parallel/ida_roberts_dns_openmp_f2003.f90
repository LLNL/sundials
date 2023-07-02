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
! This simple example problem for FIDA, due to Robertson, is from 
! chemical kinetics, and consists of the following three equations:
!
!      dy1/dt = -.04*y1 + 1.e4*y2*y3
!      dy2/dt =  .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!         0   = y1 + y2 + y3 - 1
!
! on the interval from t = 0.0 to t = 4.e10, with initial
! conditions: y1 = 1, y2 = y3 = 0.
!
! While integrating the system, we also employ the rootfinding feature
! to find the points at which y1 = 1.e-4 or at which y3 = 0.01.
!
! The problem is solved using a dense linear solver, with a
! user-supplied Jacobian. Output is printed at
! t = .4, 4, 40, ..., 4e10.
! ------------------------------------------------------------------

module dns_OMP_mod

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    integer(c_long), parameter :: neq = 3
    integer(c_long), parameter :: nout = 12

  contains

    ! ----------------------------------------------------------------
    ! fcnirob: The implicit RSH operator function
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer(c_int) function fcnres(tres, sunvec_y, sunvec_f, sunvec_r, &
            user_data) result(ierr) bind(C,name='fcnirob')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fnvector_openmp_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      real(c_double), value :: tres      ! current time
      type(N_Vector)        :: sunvec_y  ! solution N_Vector
      type(N_Vector)        :: sunvec_f  ! function N_Vector
      type(N_Vector)        :: sunvec_r  ! residual N_Vector
      type(c_ptr),    value :: user_data ! user-defined data

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yval(:)
      real(c_double), pointer :: fval(:)
      real(c_double), pointer :: rval(:)

      !======= Internals ============

      ! get data arrays from SUNDIALS vectors
      yval  => FN_VGetArrayPointer(sunvec_y)
      fval  => FN_VGetArrayPointer(sunvec_f)
      rval  => FN_VGetArrayPointer(sunvec_r)

      ! fill residual vector
      rval(1)  = -0.04d0*yval(1) + 1.0d4*yval(2)*yval(3)
      rval(2)  = -rval(1) - 3.0d7*yval(2)**2 - fval(2)
      rval(1)  = rval(1) - fval(1)
      rval(3)  = yval(1) + yval(2) + yval(3) - 1.d0

      ! return success
      ierr = 0
      return

    end function fcnres


    ! ----------------------------------------------------------------
    ! grob: The root function routine
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer(c_int) function grob(t, sunvec_y, gout, user_data) &
         result(ierr) bind(C,name='grob')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fnvector_openmp_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      real(c_double), value :: t         ! current time
      type(N_Vector)        :: sunvec_y  ! solution N_Vector
      real(c_double)        :: gout(2)   ! root function values
      type(c_ptr),    value :: user_data ! user-defined data

      ! pointers to data in SUNDIALS vectors
      real(c_double), pointer :: yval(:)

      !======= Internals ============

      ! get data array from SUNDIALS vector
      yval => FN_VGetArrayPointer(sunvec_y)

      ! fill root vector
      gout(1) = yval(1) - 0.0001d0
      gout(2) = yval(3) - 0.01d0

      ! return success
      ierr = 0
      return

    end function grob

    ! ----------------------------------------------------------------
    ! jacrob: The DAE Jacobian function
    !
    ! Return values:
    !    0 = success,
    !    1 = recoverable error,
    !   -1 = non-recoverable error
    ! ----------------------------------------------------------------
    integer(c_int) function jacrob(t, cj, sunvec_y, sunvec_f, &
         sunvec_r, sunmat_J, user_data, sunvec_t1, sunvec_t2, &
         sunvec_t3) result(ierr) bind(C,name='jacrob')

      !======= Inclusions ===========
      use, intrinsic :: iso_c_binding
      use fsundials_nvector_mod
      use fsundials_matrix_mod
      use fnvector_openmp_mod
      use fsunmatrix_dense_mod

      !======= Declarations =========
      implicit none

      ! calling variables
      real(c_double), value :: t         ! current time
      real(c_double), value :: cj        ! Jacobian scalar
      type(N_Vector)        :: sunvec_y  ! solution N_Vector
      type(N_Vector)        :: sunvec_f  ! residual N_Vector
      type(N_Vector)        :: sunvec_r  ! residual N_Vector
      type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
      type(c_ptr),    value :: user_data ! user-defined data
      type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
      type(N_Vector)        :: sunvec_t2
      type(N_Vector)        :: sunvec_t3

      ! pointers to data in SUNDIALS vector and matrix
      real(c_double), pointer :: yval(:)
      real(c_double), pointer :: J(:,:)


      !======= Internals ============

      ! get data arrays from SUNDIALS vectors
      yval => FN_VGetArrayPointer(sunvec_y)
      J(1:neq, 1:neq) => FSUNDenseMatrix_Data(sunmat_J)

      ! fill Jacobian entries
      J(1,1) = -0.04d0 - cj
      J(2,1) = 0.04d0
      J(3,1) = 1.d0
      J(1,2) = 1.d4*yval(3)
      J(2,2) = -1.d4*yval(3) - 6.d7*yval(2) - cj
      J(3,2) = 1.d0
      J(1,3) = 1.d4*yval(2)
      J(2,3) = -1.d4*yval(2)
      J(3,3) = 1.d0

      ! return success
      ierr = 0
      return

    end function jacrob

  end module dns_OMP_mod
  ! ------------------------------------------------------------------


  program main

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use, intrinsic :: omp_lib

    use fida_mod                      ! Fortran interface to IDA
    use fsundials_context_mod         ! Fortran interface to SUNContext
    use fnvector_openmp_mod           ! Fortran interface to OpenMP N_Vector
    use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
    use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
    use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
    use fsundials_matrix_mod          ! Fortran interface to generic SUNMatrix
    use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
    use fsundials_linearsolver_mod    ! Fortran interface to generic SUNLinearSolver
    use dns_OMP_mod                   ! ODE functions

    !======= Declarations =========
    implicit none

    ! local variables
    real(c_double) :: rtol, t0, tout1, tout, tret(1)
    integer(c_int) :: iout, retval, retvalr, numthreads, nrtfn, rootsfound(2)

    type(N_Vector),           pointer :: sunvec_y      ! sundials solution vector
    type(N_Vector),           pointer :: sunvec_dky    ! sundials solution vector
    type(N_Vector),           pointer :: sunvec_f      ! sundials solution vector
    type(N_Vector),           pointer :: sunvec_av     ! sundials tolerance vector
    type(SUNMatrix),          pointer :: sunmat_A      ! sundials matrix
    type(SUNLinearSolver),    pointer :: sunlinsol_LS  ! sundials linear solver
    type(c_ptr)                       :: ida_mem    ! IDA memory
    type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

    ! solution and tolerance vectors, neq is set in the dae_mod module
    real(c_double) :: yval(neq), fval(neq), avtol(neq), dkyval(neq)

    !======= Internals ============

    retval = FSUNContext_Create(c_null_ptr, sunctx)

    ! initialize solution vectors and tolerances
    yval(1) = 1.d0
    yval(2) = 0.d0
    yval(3) = 0.d0

    fval(1) = -0.04d0
    fval(2) = 0.04d0
    fval(3) = 0.d0

    rtol = 1.d-4

    avtol(1) = 1.d-6
    avtol(2) = 1.d-10
    avtol(3) = 1.d-6

    numthreads = OMP_get_num_threads()

    ! create serial vectors
    sunvec_y => FN_VMake_OpenMP(neq, yval, numthreads, sunctx)
    if (.not. associated(sunvec_y)) then
       print *, 'ERROR: sunvec = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    sunvec_f => FN_VMake_OpenMP(neq, fval, numthreads, sunctx)
    if (.not. associated(sunvec_f)) then
       print *, 'ERROR: sunvec = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    sunvec_av => FN_VMake_OpenMP(neq, avtol, numthreads, sunctx)
    if (.not. associated(sunvec_av)) then
       print *, 'ERROR: sunvec = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! set integration limits
    t0 = 0.d0
    tout1 = 0.4d0

    call PrintHeader(rtol, avtol, yval)

    ! Call FIDACreate and FIDAInit to create and initialize IDA memory
    ida_mem = FIDACreate(sunctx)
    if (.not. c_associated(ida_mem)) print *, 'ERROR: ida_mem = NULL'

    ! parallel? shared(rtol, avtol)

    retval = FIDAInit(ida_mem, c_funloc(fcnres), t0, sunvec_y, sunvec_f)
    if (retval /= 0) then
       print *, 'Error in FIDAInit, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! Call FIDASVtolerances to set tolerances
    retval = FIDASVtolerances(ida_mem, rtol, sunvec_av)
    if (retval /= 0) then
       print *, 'Error in FIDASVtolerances, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! Call FIDARootInit to specify the root function grob with 2 components
    nrtfn = 2
    retval = FIDARootInit(ida_mem, nrtfn, c_funloc(grob))
    if (retval /= 0) then
       print *, 'Error in FIDARootInit, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! end parallel?

    ! Create dense SUNMatrix for use in linear solves
    sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
    if (.not. associated(sunmat_A)) then
       print *, 'ERROR: sunmat = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! Create dense SUNLinearSolver object
    sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
    if (.not. associated(sunlinsol_LS)) then
       print *, 'ERROR: sunlinsol = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! Attach the matrix and linear solver
    retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
    if (retval /= 0) then
       print *, 'Error in FIDASetLinearSolver, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! parallel?

    ! Set the user-supplied Jacobian routine
    retval = FIDASetJacFn(ida_mem, c_funloc(jacrob))
    if (retval /= 0) then
       print *, 'Error in FIDASetJacFn, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! end parallel?

    ! In loop, call IDASolve, print results, and test for error.

    iout = 0
    tout = tout1
    do while(iout < nout)

       ! parallel?
       retval = FIDASolve(ida_mem, tout, tret(1), sunvec_y, sunvec_f, IDA_NORMAL)
       if (retval < 0) then
          print *, 'Error in FIDASolve, retval = ', retval, '; halting'
          call FIDAFree(ida_mem)
          stop 1
       end if
       ! end parallel?

       call PrintOutput(ida_mem, tret(1), yval)

       if (retval .eq. IDA_ROOT_RETURN) then
          ! parallel?
          retvalr = FIDAGetRootInfo(ida_mem, rootsfound)
          if (retvalr < 0) then
             print *, 'Error in FIDAGetRootInfo, retval = ', retval, '; halting'
             call FIDAFree(ida_mem)
             stop 1
          end if
          ! end parallel?
          print '(a,2(i2,2x))', "    rootsfound[] = ", rootsfound(1), rootsfound(2)
       end if

       if (retval .eq. IDA_SUCCESS) then
          iout = iout + 1
          tout = tout * 10.d0
       end if
    end do

    sunvec_dky => FN_VMake_OpenMP(neq, dkyval, numthreads, sunctx)
    if (.not. associated(sunvec_dky)) then
       print *, 'ERROR: sunvec = NULL'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! parallel?

    ! find and print derivative at tret(1)
    retval = FIDAGetDky(ida_mem, tret(1), 1, sunvec_dky)
    if (retval /= 0) then
       print *, 'Error in IDAGetDky'
       call FIDAFree(ida_mem)
       stop 1
    end if

    ! end parallel?

    print *, " "
    print *, "------------------------------------------------------"
    print *, "     Final      y1'          y2'          y3'"
    print *, "------------------------------------------------------"
    print '(13x,3(es12.4,1x))', dkyval

    call PrintFinalStats(ida_mem)

    ! free memory
    call FIDAFree(ida_mem)
    retval = FSUNLinSolFree(sunlinsol_LS)
    call FSUNMatDestroy(sunmat_A)
    call FN_VDestroy_OpenMP(sunvec_y)
    call FN_VDestroy_OpenMP(sunvec_dky)
    call FN_VDestroy_OpenMP(sunvec_av)
    retval = FSUNContext_Free(sunctx)

  end program main


  ! ----------------------------------------------------------------
  ! PrintHeader: prints first lines of output (problem description)
  ! ----------------------------------------------------------------
  subroutine PrintHeader(rtol, avtol, y)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use dns_OMP_mod

    !======= Declarations =========
    implicit none

    ! calling variable
    real(c_double) :: rtol
    real(c_double) :: avtol(neq)
    real(c_double) :: y(neq)

    !======= Internals ============

    print *, " "
    print *, "ida_roberts_dns_openmp_f2003.f90: Robertson IDA ODE serial example problem for IDA"
    print *, "         Three equation chemical kinetics problem."
    print *, " "
    print *, "Linear solver: DENSE, with user-supplied Jacobian."
    print '(a,f6.4,a,3(es7.0,1x))', "Tolerance parameters:  rtol = ",rtol,"   atol = ", avtol
    print '(a,3(f5.2,1x),a)', "Initial conditions y0 = (",y,")"
    print *, "Constraints and id not used."
    print *, " "
    print *, "------------------------------------------------------------------------"
    print *, "   t            y1           y2           y3       | nst   k     h"
    print *, "------------------------------------------------------------------------"

    return
  end subroutine PrintHeader


  ! ----------------------------------------------------------------
  ! PrintOutput
  ! ----------------------------------------------------------------
  subroutine PrintOutput(ida_mem, t, y)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fida_mod
    use dns_OMP_mod

    !======= Declarations =========
    implicit none

    ! calling variable
    type(c_ptr)    :: ida_mem
    real(c_double) :: t, y(neq)

    ! internal variables
    integer(c_int)  :: retval, kused(1)
    integer(c_long) :: nst(1)
    real(c_double)  :: hused(1)

    !======= Internals ============

    retval = FIDAGetNumSteps(ida_mem, nst)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetLastOrder(ida_mem, kused)
    if (retval /= 0) then
       print *, 'Error in FIDAGetLastOrder, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
     end if

    retval = FIDAGetLastStep(ida_mem, hused)
    if (retval /= 0) then
       print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    print '(es12.4,1x,3(es12.4,1x),a,2(i3,2x),es12.4)', &
         t, y(1), y(2), y(3), "| ", nst, kused(1), hused(1)

  end subroutine PrintOutput


  ! ----------------------------------------------------------------
  ! PrintFinalStats
  !
  ! Print IDASOL statstics to standard out
  ! ----------------------------------------------------------------
  subroutine PrintFinalStats(ida_mem)

    !======= Inclusions ===========
    use iso_c_binding
    use fida_mod

    !======= Declarations =========
    implicit none

    type(c_ptr), intent(inout) :: ida_mem ! solver memory structure

    integer(c_int)  :: retval          ! error flag

    integer(c_long) :: nsteps(1)     ! num steps
    integer(c_long) :: nre(1)        ! num residual function evals
    integer(c_long) :: netfails(1)   ! num error test fails

    integer(c_long) :: nniters(1)    ! nonlinear solver iterations
    integer(c_long) :: nncfails(1)   ! nonlinear solver fails
    integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations
    integer(c_long) :: ngevals(1)  ! number of root function evaluations

    !======= Internals ============

    retval = FIDAGetNumSteps(ida_mem, nsteps)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumResEvals(ida_mem, nre)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumRhsEvals, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumErrTestFails(ida_mem, netfails)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumErrTestFails, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumNonlinSolvIters(ida_mem, nniters)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumNonlinSolvIters, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumNonlinSolvConvFails(ida_mem, nncfails)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumGEvals(ida_mem, ngevals)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumJacEvals, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    retval = FIDAGetNumJacEvals(ida_mem, njacevals)
    if (retval /= 0) then
       print *, 'Error in FIDAGetNumJacEvals, retval = ', retval, '; halting'
       call FIDAFree(ida_mem)
       stop 1
    end if

    print *, ' '
    print *, ' General Solver Stats:'
    print '(4x,A,i9)'    ,'Total internal steps taken    =',nsteps
    print '(4x,A,i9)'    ,'Total residual function calls =',nre
    print '(4x,A,i9)'    ,'Total Jacobian function calls =',njacevals
    print '(4x,A,i9)'    ,'Total root fcn function calls =',ngevals
    print '(4x,A,i9)'    ,'Num error test failures       =',netfails
    print '(4x,A,i9)'    ,'Num nonlinear solver iters    =',nniters
    print '(4x,A,i9)'    ,'Num nonlinear solver fails    =',nncfails
    print *, ' '

    return

  end subroutine PrintFinalStats
