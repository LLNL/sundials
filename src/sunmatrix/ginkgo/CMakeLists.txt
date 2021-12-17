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

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_GINKGO with ${SUNDIALS_GINKGO_BACKENDS} backend(s)\n\")")

# if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
#   set_source_files_properties(sunmatrix_ginkgo.cpp PROPERTIES LANGUAGE CUDA)
#   set(_libs_needed sundials_nveccuda ${CUDA_CUBLAS_LIBRARIES})
# elseif(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
#   set_source_files_properties(sunmatrix_ginkgo.cpp PROPERTIES LANGUAGE CXX)
#   set(_libs_needed sundials_nvechip hip::device)
# endif()

# Add the sunmatrix_ginkgo library
sundials_add_library(sundials_sunmatrixginkgo
  SOURCES
    sunmatrix_ginkgodense.cpp
    sunmatrix_ginkgocsr.cpp
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgo.h
  INCLUDE_SUBDIR
    sunmatrix
  OBJECT_LIBRARIES
    sundials_generic_obj
  # INCLUDE_DIRECTORIES
    # PUBLIC ${MAGMA_INCLUDE_DIR}
  LINK_LIBRARIES
    # PUBLIC SUNDIALS::GINKGO ${_libs_needed}
    PUBLIC Ginkgo::ginkgo
  OUTPUT_NAME
    sundials_sunmatrixginkgo
  VERSION
    ${sunmatrixlib_VERSION}
  SOVERSION
    ${sunmatrixlib_SOVERSION}
)

message(STATUS "Added SUNMATRIX_GINKGO module with ${SUNDIALS_GINKGO_BACKENDS} backend(s)")
