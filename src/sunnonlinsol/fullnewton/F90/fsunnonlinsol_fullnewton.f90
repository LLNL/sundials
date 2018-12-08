! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2017, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! -----------------------------------------------------------------
! This file contains a Fortran module for interfacing directly with
! the SUNDIALS Full Newton iteration nonlinear solver using the
! ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunnonlinsol_fullnewton_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNNonlinSol_FullNewton(y) &
         bind(C,name='SUNNonlinSol_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
     end function FSUNNonlinSol_FullNewton

     ! =================================================================
     ! Destructors
     ! =================================================================

     integer(c_int) function FSUNNonlinSolFree_FullNewton(NLS) &
         bind(C,name='SUNNonlinSolFree_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolFree_FullNewton

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNNonlinSolGetType_FullNewton(NLS) &
         bind(C,name='SUNNonlinSolGetType_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolGetType_FullNewton

     integer(c_int) function FSUNNonlinSolInitialize_FullNewton(NLS) &
         bind(C,name='SUNNonlinSolInitialize_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolInitialize_FullNewton

     integer(c_int) function FSUNNonlinSolSolve_FullNewton(NLS, y0, y, w, tol, &
                                                           callSetup, mem) &
         bind(C,name='SUNNonlinSolSolve_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_ptr),    value :: y0
       type(c_ptr),    value :: y
       type(c_ptr),    value :: w
       real(c_double), value :: tol
       integer(c_int), value :: callSetup
       type(c_ptr),    value :: mem
     end function FSUNNonlinSolSolve_FullNewton
    
    ! =================================================================
    ! Set functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolSetSysFn_FullNewton(NLS, SysFn) &
         bind(C,name='SUNNonlinSolSetSysFn_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: SysFn
     end function FSUNNonlinSolSetSysFn_FullNewton

     integer(c_int) function FSUNNonlinSolSetLSetupFn_FullNewton(NLS, LSetupFn) &
         bind(C,name='SUNNonlinSolSetLSetupFn_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: LSetupFn
     end function FSUNNonlinSolSetLSetupFn_FullNewton

     integer(c_int) function FSUNNonlinSolSetLSolveFn_FullNewton(NLS, LSolveFn) &
         bind(C,name='SUNNonlinSolSetLSolveFn_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: LSolveFn
     end function FSUNNonlinSolSetLSolveFn_FullNewton

     integer(c_int) function FSUNNonlinSolSetConvTestFn_FullNewton(NLS, CTestFN) &
         bind(C,name='SUNNonlinSolSetConvTestFn_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: CTestFN
     end function FSUNNonlinSolSetConvTestFn_FullNewton

     integer(c_int) function FSUNNonlinSolSetMaxIters_FullNewton(NLS, maxiters) &
         bind(C,name='SUNNonlinSolSetMaxIters_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_int), value :: maxiters
     end function FSUNNonlinSolSetMaxIters_FullNewton
 
    ! =================================================================
    ! Get functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolGetNumIters_FullNewton(NLS, niters) &
         bind(C,name='SUNNonlinSolGetNumIters_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_long)       :: niters
     end function FSUNNonlinSolGetNumIters_FullNewton

     integer(c_int) function FSUNNonlinSolGetCurIter_FullNewton(NLS, iter) &
         bind(C,name='SUNNonlinSolGetCurIter_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       integer(c_int)       :: iter
     end function FSUNNonlinSolGetCurIter_FullNewton
     
     integer(c_int) function FSUNNonlinSolGetSysFn_FullNewton(NLS, SysFn) &
         bind(C,name='SUNNonlinSolGetSysFn_FullNewton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       type(c_funptr)       :: SysFn
     end function FSUNNonlinSolGetSysFn_FullNewton

   end interface

end module fsunnonlinsol_fullnewton_mod
