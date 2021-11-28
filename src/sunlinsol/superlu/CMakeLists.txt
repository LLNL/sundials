# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the superlu SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLU\n\")")

# Add the library
sundials_add_library(sundials_sunlinsolsuperlu
  SOURCES
    sunlinsol_superlu.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_superlu.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC
    sundials_sunmatrixsparse
    SUNDIALS::SUPERLU
  OUTPUT_NAME
    sundials_sunlinsolsuperlu
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
    ${sunlinsollib_VERSION}
)

message(STATUS "Added SUNLINSOL_SUPERLU module")

# If FCMIX is enabled, build and install the Fortran77 library
if(BUILD_FORTRAN77_INTERFACE)
  sundials_add_library(sundials_fsunlinsolsuperlu
    SOURCES
      fsunlinsol_superlu.c
    LINK_LIBRARIES
      PUBLIC sundials_sunlinsolsuperlu
      PRIVATE sundials_fnvecserial sundials_fsunmatrixsparse
    OUTPUT_NAME
      sundials_fsunlinsolsuperlu
    VERSION
      ${sunlinsollib_VERSION}
    SOVERSION
      ${sunlinsollib_SOVERSION}
  )
  message(STATUS "Added SUNLINSOL_SUPERLU F77 interface")
endif()
