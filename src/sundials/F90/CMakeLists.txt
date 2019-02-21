# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the F2003 SUNDIALS object library
# ---------------------------------------------------------------

set(sundials_SOURCES
  fnvector_mod.f90
  fsundials_types.f90
  fsunmatrix_mod.f90
  fsunlinsol_mod.f90
  fsunnonlinsol_mod.f90)

if(BUILD_STATIC_LIBS)
  add_library(sundials_fgeneric_static_obj OBJECT ${sundials_SOURCES})
  set_target_properties(sundials_fgeneric_static_obj PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_STATIC)
endif()

if(BUILD_SHARED_LIBS)
  add_library(sundials_fgeneric_shared_obj OBJECT ${sundials_SOURCES})
  set_target_properties(sundials_fgeneric_shared_obj PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}_SHARED)
endif()

