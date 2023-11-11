# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the superlu-dist SUNLinearSolver library
# -----------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SUPERLUDIST\n\")")

# Add the library
sundials_add_library(sundials_sunlinsolsuperludist
  SOURCES
    sunlinsol_superludist.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_superludist.h
  INCLUDE_SUBDIR
    sunlinsol
  LINK_LIBRARIES
    PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES
    PUBLIC
    sundials_sunmatrixslunrloc
    SUNDIALS::SUPERLUDIST
    $<IF:$<BOOL:${SUPERLUDIST_OpenMP}>,OpenMP::OpenMP_C,>
    MPI::MPI_C
  COMPILE_OPTIONS
    PRIVATE ${_compile_options}
  OUTPUT_NAME
    sundials_sunlinsolsuperludist
  VERSION
    ${sunlinsollib_VERSION}
  SOVERSION
  ${sunlinsollib_SOVERSION}
)

message(STATUS "Added SUNLINSOL_SUPERLUDIST module")
