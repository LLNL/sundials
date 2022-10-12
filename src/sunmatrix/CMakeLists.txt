# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# matrix level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required native matrices
add_subdirectory(band)
add_subdirectory(dense)
add_subdirectory(sparse)

# optional TPL matrices
if(BUILD_SUNMATRIX_CUSPARSE)
  add_subdirectory(cusparse)
endif()

if(BUILD_SUNMATRIX_GINKGO)
  install(CODE "MESSAGE(\"\nInstall SUNMATRIX_GINKGO\n\")")
  message(STATUS "Added BUILD_SUNMATRIX_GINKGO module")
  add_library(sundials_sunmatrixginko INTERFACE)
  target_include_directories(sundials_sunmatrixginko INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${SUNDIALS_SOURCE_DIR}/src/sundials>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  install(FILES ${PROJECT_SOURCE_DIR}/include/sunmatrix/sunmatrix_ginkgo.hpp
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sunmatrix")
  set(_SUNDIALS_INSTALLED_COMPONENTS "sunmatrixginkgo;${_SUNDIALS_INSTALLED_COMPONENTS}" CACHE INTERNAL "" FORCE)
endif()

if(BUILD_SUNMATRIX_KOKKOSDENSE)
  add_subdirectory(kokkos)
endif()

if(BUILD_SUNMATRIX_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNMATRIX_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNMATRIX_SLUNRLOC)
  add_subdirectory(slunrloc)
endif()
