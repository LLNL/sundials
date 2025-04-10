# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda cuSolverSp SUNLinearSolver
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_CUSOLVERSP\n\")")

sundials_add_library(
  sundials_sunlinsolcusolversp
  SOURCES sunlinsol_cusolversp_batchqr.cu
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_cusolversp_batchqr.h
  INCLUDE_SUBDIR sunlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES PUBLIC sundials_sunmatrixcusparse CUDA::cusolver PRIVATE
                 CUDA::cusparse
  OUTPUT_NAME sundials_sunlinsolcusolversp
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_SOVERSION})

message(STATUS "Added SUNLINSOL_CUSOLVERSP module")
