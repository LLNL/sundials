# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the LapackDense SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKDENSE\n\")")

# Add the library
sundials_add_library(sundials_sunlinsollapackdense
  SOURCES
    sunlinsol_lapackdense.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackdense.h
  INCLUDE_SUBDIR
    sunlinsol
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC sundials_sunmatrixdense "${LAPACK_LIBRARIES}"
  OUTPUT_NAME
    sundials_sunlinsollapackdense
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
  ${sunlinsollib_SOVERSION}
)

message(STATUS "Added SUNLINSOL_LAPACKDENSE module")

# Add module if F2003 interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
