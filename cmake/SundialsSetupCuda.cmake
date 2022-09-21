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

set(DOCSTR "The CUDA standard to use if CUDA is enabled (11, 14, 17, 20)")
sundials_option(CMAKE_CUDA_STANDARD STRING "${DOCSTR}" "${CMAKE_CXX_STANDARD}"
                OPTIONS "11;14;17;20")
message(STATUS "CUDA standard set to ${CMAKE_CUDA_STANDARD}")

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

# ===============================================================
# Enable CUDA lang and find the CUDA libraries.
# ===============================================================

enable_language(CUDA)
set(CUDA_FOUND TRUE)

# Need this as long as CUDA libraries like cuSOLVER are not available
# through some other way.
find_package(CUDA REQUIRED)

# Hide legacy FindCUDA variables
get_cmake_property(_variables VARIABLES)
foreach(_var ${_variables})
  if("${_var}" MATCHES "^CUDA_[A-z]+_LIBRARY")
    # do nothing
  elseif("${_var}" MATCHES "^CUDA_.*")
    mark_as_advanced(${_var})
  endif()
endforeach()

# Make the CUDA_rt_LIBRARY advanced like the other CUDA_*_LIBRARY variables
mark_as_advanced(FORCE CUDA_rt_LIBRARY)

# Show CUDA flags
mark_as_advanced(CLEAR CMAKE_CUDA_FLAGS)

# ===============================================================
# Print out information about CUDA.
# ===============================================================

message(STATUS "CUDA Version:               ${CUDA_VERSION_STRING}")
message(STATUS "CUDA Architectures:         ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CUDA Compiler:              ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA Host Compiler:         ${CMAKE_CUDA_HOST_COMPILER}")
message(STATUS "CUDA Include Path:          ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Libraries:             ${CUDA_LIBRARIES}")
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
