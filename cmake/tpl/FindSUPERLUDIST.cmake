# ---------------------------------------------------------------
# Programmer(s): Cody Balos @ LLNL
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
# SuperLUDIST find module that creates an imported target for
# SuperLU_DIST. The target is SUNDIALS::SUPERLUDIST.
#
# The module has two modes: one 'automatic' mode that uses pkg-conf
# to locate and determine the SuperLU DIST libraries and includes,
# and one 'manual' mode which requires the libraries and includes
# to be specified.
#
# If SUPERLUDIST_LIBRARIES is set, then the 'manual' mode is chosen
# and the SUPERLUDIST_INCLUDE_DIR variable must also be set.
#
# The variable SUPERLUDIST_DIR or SUPERLUDIST_LIBRARY_DIR can be
# used to control where the module looks for the library in
# 'automatic' mode.
#
# This module also defines additional variables:
#
#   SUPERLUDIST_FOUND      - the SuperLU_DIST libraries were found
#   SUPERLUDIST_INDEX_SIZE - the bit width of indices in SUPERLUDIST
#   SUPERLUDIST_CUDA       - SuperLU_DIST has CUDA support
#   SUPERLUDIST_ROCM       - SuperLU_DIST has ROCm support
# ---------------------------------------------------------------

# Check if SUPERLUDIST_LIBRARIES contains the superlu_dist
# library as well as TPLs. If so, extract it into the
# SUPERLUDIST_LIBRARY variable.
if(NOT SUPERLUDIST_LIBRARIES)
  find_package(PkgConfig REQUIRED)
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${SUPERLUDIST_DIR}")
  pkg_search_module(SUPERLUDIST REQUIRED superlu_dist>=6.1.1)
  pkg_get_variable(SUPERLUDIST_INCLUDE_DIR superlu_dist includedir)

  # find the library configuration file
  set(SUPERLUDIST_CUDA FALSE CACHE BOOL "SuperLU DIST was built with CUDA support")
  set(SUPERLUDIST_ROCM FALSE CACHE BOOL "SuperLU DIST was built with ROCm support")
  if(SUPERLUDIST_INCLUDE_DIR)
    find_file(SUPERLUDIST_CONFIG_PATH superlu_dist_config.h PATHS ${SUPERLUDIST_INCLUDE_DIR})
    mark_as_advanced(FORCE SUPERLUDIST_CONFIG_PATH)
    if(SUPERLUDIST_VERSION VERSION_GREATER_EQUAL "8.0.0")
        file(STRINGS ${SUPERLUDIST_CONFIG_PATH} _index_size_64 REGEX "#define XSDK_INDEX_SIZE 64")
        file(STRINGS ${SUPERLUDIST_CONFIG_PATH} _index_size_32 REGEX "#undef XSDK_INDEX_SIZE")
        if(_index_size_64)
          set(SUPERLUDIST_INDEX_SIZE 64 CACHE STRING "SuperLU DIST index size (bit width)" FORCE)
        else()
          set(SUPERLUDIST_INDEX_SIZE 32 CACHE STRING "SuperLU DIST index size (bit width)" FORCE)
        endif()
        mark_as_advanced(FORCE SUPERLUDIST_INDEX_SIZE)
    endif()
    file(STRINGS ${SUPERLUDIST_CONFIG_PATH} _strings_have_cuda REGEX "HAVE_CUDA")
    string(REGEX MATCH "TRUE|FALSE" _has_cuda "${_strings_have_cuda}")
    file(STRINGS ${SUPERLUDIST_CONFIG_PATH} _strings_have_rocm REGEX "HAVE_ROCM")
    string(REGEX MATCH "TRUE|FALSE" _has_rocm "${_strings_have_rocm}")
    set(SUPERLUDIST_CUDA ${_has_cuda} CACHE BOOL "SuperLU DIST was built with CUDA support" FORCE)
    set(SUPERLUDIST_ROCM ${_has_rocm} CACHE BOOL "SuperLU DIST was built with ROCm support" FORCE)
    unset(_has_cuda)
    unset(_has_rocm)
    mark_as_advanced(SUPERLUDIST_CUDA)
    mark_as_advanced(SUPERLUDIST_ROCM)
  endif()

  set(SUPERLUDIST_INTERFACE_LINK_LIBRARIES "${SUPERLUDIST_LINK_LIBRARIES}" CACHE INTERNAL "")

  # add libraries for OpenMP support
  if(SUPERLUDIST_OpenMP)
    list(APPEND SUPERLUDIST_INTERFACE_LINK_LIBRARIES OpenMP::OpenMP_CXX)
  endif()
  # add libraries for GPU support
  if(SUPERLUDIST_CUDA)
    list(APPEND SUPERLUDIST_INTERFACE_LINK_LIBRARIES CUDA::cudart CUDA::cublas)
  endif()
  if(SUPERLUDIST_HIP)
     list(APPEND SUPERLUDIST_INTERFACE_LINK_LIBRARIES roc::hipblas)
  endif()
endif()

set(SUPERLUDIST_INTERFACE_LINK_LIBRARIES "${SUPERLUDIST_INTERFACE_LINK_LIBRARIES}" CACHE INTERNAL "")

# set package variables including SUPERLUDIST_FOUND
find_package_handle_standard_args(SUPERLUDIST
  REQUIRED_VARS
    SUPERLUDIST_INTERFACE_LINK_LIBRARIES
    SUPERLUDIST_INCLUDE_DIRS
    SUPERLUDIST_INDEX_SIZE
  VERSION_VAR
    SUPERLUDIST_VERSION
  )

# Create target for SuperLU_DIST
if(SUPERLUDIST_FOUND)

  if(NOT TARGET SUNDIALS::SUPERLUDIST)
    add_library(SUNDIALS::SUPERLUDIST INTERFACE IMPORTED)
  endif()

  set_target_properties(SUNDIALS::SUPERLUDIST PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUDIST_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${SUPERLUDIST_INTERFACE_LINK_LIBRARIES}")

endif()
