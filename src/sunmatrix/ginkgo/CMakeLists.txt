# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_GINKGO with ${GINKGO_BACKENDS} backend(s)\n\")")

# Add the sunmatrix_ginkgo library
sundials_add_library(sundials_sunmatrixginkgo
  SOURCES
    sunmatrix_ginkgodense.cpp
    sunmatrix_ginkgocsr.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgo.hpp
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgoblock.hpp
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
  LINK_LIBRARIES
    PUBLIC Ginkgo::ginkgo
  OUTPUT_NAME
    sundials_sunmatrixginkgo
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_GINKGO module with ${GINKGO_BACKENDS} backend(s)")
