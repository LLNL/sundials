# ---------------------------------------------------------------
# Programmer(s): Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the parhyp NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_PARHYP with ${SUNDIALS_HYPRE_BACKENDS} backend\n\")")

if(SUNDIALS_HYPRE_BACKENDS MATCHES "SERIAL")
  # Do nothing; Don't use sumemoryhelper, don't change language
elseif(SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
  set_source_files_properties(nvector_parhyp.c PROPERTIES LANGUAGE CUDA)
  set(_sunmemlib sundials_sunmemcuda_obj) # Use CUDA sunmemoryhelper
  set(_hypre_using_cuda TRUE)
#   set(_include_vectorkernels_cuh TRUE)    # Include src/nvector/cuda/VectorKernels.cuh
  set(_cxx_std cxx_std_14)
elseif(SUNDIALS_HYPRE_BACKENDS MATCHES "HIP")
  set_source_files_properties(nvector_parhyp.c PROPERTIES LANGUAGE HIP)
  set(_sunmemlib sundials_sunmemhip_obj)  # Use HIP sunmemoryhelper
  set(_hypre_using_hip TRUE)
#   set(_include_vectorkernels_hip TRUE)    # Include src/nvector/hip/VectorKernels.hip.hpp
  set(_hip_lib_or_not hip::device)        # Require hip::device library
  set(_cxx_std cxx_std_14)
endif()

# Create the library
sundials_add_library(sundials_nvecparhyp
  SOURCES
    nvector_parhyp.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_parhyp.h
  INCLUDE_SUBDIR
    nvector
  OBJECT_LIBRARIES
    sundials_generic_obj
    ${_sunmemlib}
  LINK_LIBRARIES
    PUBLIC
    MPI::MPI_C
    SUNDIALS::HYPRE
    $<$<BOOL:_hypre_using_cuda>:${std::cstdio}>
    $<$<BOOL:_hypre_using_cuda>:${CMAKE_CPP_STANDARD_LIBRARIES}>
    $<$<BOOL:_hypre_using_cuda>:${CMAKE_CUDA_STANDARD_LIBRARIES}>
    ${_hip_lib_or_not}
  INCLUDE_DIRECTORIES
    PUBLIC
    $<$<BOOL:_hypre_using_cuda>:${SUNDIALS_SOURCE_DIR}/src/nvector/cuda>
    $<$<BOOL:_hypre_using_cuda>:${CMAKE_CPP_STANDARD_INCLUDE_DIRECTORIES}>
    $<$<BOOL:_hypre_using_cuda>:${CMAKE_CUDA_STANDARD_INCLUDE_DIRECTORIES}>
    $<$<BOOL:_hypre_using_hip>:${SUNDIALS_SOURCE_DIR}/src/nvector/hip>
  COMPILE_DEFINITIONS
    PUBLIC
    USE_CUDA
#TODO Is this needed? -jsdomine
#  COMPILE_OPTIONS
#    PUBLIC
#    $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
  COMPILE_FEATURES
    PUBLIC
    ${_cxx_std}
  PROPERTIES
    $<$<BOOL:_hypre_using_cuda>:LANGUAGE CUDA>
  OUTPUT_NAME
    sundials_nvecparhyp
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_PARHYP module with ${SUNDIALS_HYPRE_BACKENDS} backend")
