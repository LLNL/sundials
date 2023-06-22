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

if(SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
  set_source_files_properties(nvector_parhyp.c PROPERTIES LANGUAGE CUDA)
  set(_memhelper sundials_sunmemcuda_obj)
  set(_cxx_std cxx_std_14)
  set(_extra_includes ${CMAKE_SOURCE_DIR}/src/nvector/cuda)
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
    ${_memhelper}
  LINK_LIBRARIES
    PUBLIC
    MPI::MPI_C
    SUNDIALS::HYPRE
  INCLUDE_DIRECTORIES
    PRIVATE
    ${_extra_includes}
  COMPILE_FEATURES
    PUBLIC
    ${_cxx_std}
  OUTPUT_NAME
    sundials_nvecparhyp
  VERSION
    ${nveclib_VERSION}
  SOVERSION
    ${nveclib_SOVERSION}
)

message(STATUS "Added NVECTOR_PARHYP module with ${SUNDIALS_HYPRE_BACKENDS} backend")
