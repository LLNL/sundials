# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the Arnoldi Iteration
# SUNDomEigEstimator library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNDOMEIGEST_ARNI\n\")")

# Workaround bug in older versions of CMake where the BLAS::BLAS target, which
# LAPACK::LAPACK depends on, is not picked up by LINK_LIBRARIES
set(_lapack_targets LAPACK::LAPACK)
if(CMAKE_VERSION VERSION_LESS 3.20)
  list(APPEND _lapack_targets BLAS::BLAS)
endif()

# Add the sundomeigest_arni library
sundials_add_library(
  sundials_sundomeigestarni
  SOURCES sundomeigest_arni.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sundomeigest/sundomeigest_arni.h
  INCLUDE_SUBDIR sundomeigest
  LINK_LIBRARIES PUBLIC sundials_core ${_lapack_targets}
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sundomeigestarni
  VERSION ${sundomeigestlib_VERSION}
  SOVERSION ${sundomeigestlib_SOVERSION})

message(STATUS "Added SUNDOMEIGEST_ARNI module")
