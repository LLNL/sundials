# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Setup the CUDA languge and CUDA libraries.
# ---------------------------------------------------------------

# For CUDA support, require CMake 3.18 so we can use FindCUDAToolkit
# FindCUDAToolkit was introduced in 3.17, but 3.18 fixes a lot
# of issues with it and CUDA as a native language.
cmake_minimum_required(VERSION 3.18.0)

# ===============================================================
# Configure options needed prior to enabling the CUDA language
# ===============================================================

if(NOT CMAKE_CUDA_HOST_COMPILER)
  # If a user did not provide the host compiler, then we
  # assume that they want to use the CXX compiler that was set.
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "NVCC host compiler")
endif()

# ===============================================================
# Configure the CUDA flags
# ===============================================================

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")

if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
  if(CMAKE_CUDA_ARCHITECTURES)
    foreach(arch ${CMAKE_CUDA_ARCHITECTURES})
      # Remove real/virtual specifiers
      string(REGEX MATCH "[0-9]+" arch_name "${arch}")
      string(APPEND _nvcc_arch_flags " -gencode=arch=compute_${arch_name},code=sm_${arch_name}")
    endforeach()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${_nvcc_arch_flags}")
  endif()
endif()

if( (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    OR (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    AND (CMAKE_SYSTEM_PROCESSOR MATCHES ppc64le) )
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-mno-float128 _hasflag)
  if(_hasflag)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-mno-float128")
  endif()
endif()

# Need c++14 for the CUDA compiler check.
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++14")

# ===============================================================
# Enable CUDA lang and find the CUDA libraries.
# ===============================================================

enable_language(CUDA)
set(CUDA_FOUND TRUE)

find_package(CUDAToolkit REQUIRED)

# Show CUDA flags
mark_as_advanced(CLEAR CMAKE_CUDA_FLAGS)

# We need c++14 for the CUDA compiler check, but if we don't remove it,
# then we will get a redefinition error. CMAKE_CUDA_STANDARD ends up
# setting the proper version.
if(CMAKE_CUDA_FLAGS)
  STRING(REPLACE "-std=c++14" " " CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS})
endif()
set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})

# ===============================================================
# Print out information about CUDA.
# ===============================================================

message(STATUS "CUDA Toolkit Version:       ${CUDAToolkit_VERSION}")
message(STATUS "CUDA Architectures:         ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CUDA Compiler:              ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA Host Compiler:         ${CMAKE_CUDA_HOST_COMPILER}")
message(STATUS "CUDA Toolkit Includes:      ${CUDAToolkit_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directory:     ${CUDAToolkit_LIBRARY_DIR}")
message(STATUS "CUDA Compile Flags:         ${CMAKE_CUDA_FLAGS}")
message(STATUS "CUDA Link Flags:            ${CMAKE_CUDA_LINK_FLAGS}")
message(STATUS "CUDA Link Executable:       ${CMAKE_CUDA_LINK_EXECUTABLE}")
message(STATUS "CUDA Separable Compilation: ${CMAKE_CUDA_SEPARABLE_COMPILATION}")


# ===============================================================
# Configure compiler for installed examples
# ===============================================================

if((SUNDIALS_BUILD_WITH_PROFILING OR SUNDIALS_LOGGING_ENABLE_MPI) AND ENABLE_MPI)
  set(_EXAMPLES_CUDA_HOST_COMPILER "${MPI_CXX_COMPILER}" CACHE INTERNAL "${lang} compiler for installed examples")
else()
  set(_EXAMPLES_CUDA_HOST_COMPILER "${CMAKE_CUDA_HOST_COMPILER}" CACHE INTERNAL "${lang} compiler for installed examples")
endif()