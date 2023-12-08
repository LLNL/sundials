# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU,
#                and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the band SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_BAND\n\")")

# Add the library
sundials_add_library(sundials_sunlinsolband
  SOURCES
    sunlinsol_band.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_band.h
  INCLUDE_SUBDIR
    sunlinsol
  LINK_LIBRARIES
    PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES
    PUBLIC sundials_sunmatrixband
  OUTPUT_NAME
    sundials_sunlinsolband
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
  ${sunlinsollib_SOVERSION}
)

message(STATUS "Added SUNLINSOL_BAND module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(fmod)
endif()
