# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the CUDA cuSPARSE SUNMatrix
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_CUSPARSE\n\")")

# Add the library
sundials_add_library(sundials_sunmatrixcusparse
  SOURCES
    sunmatrix_cusparse.cu
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_cusparse.h
  INCLUDE_SUBDIR
    sunmatrix
  LINK_LIBRARIES
    PUBLIC sundials_core
  OBJECT_LIBRARIES
    sundials_sunmemcuda_obj
  LINK_LIBRARIES
    PUBLIC CUDA::cusparse CUDA::cusolver
  OUTPUT_NAME
    sundials_sunmatrixcusparse
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_CUSPARSE module")
