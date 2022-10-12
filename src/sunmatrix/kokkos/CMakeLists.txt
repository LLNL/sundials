# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the Kokkos SUNMatrix
# ------------------------------------------------------------------------------

if(BUILD_SUNMATRIX_KOKKOSDENSE)
  install(CODE "MESSAGE(\"\nInstall SUNMATRIX_KOKKOSDENSE\n\")")

  sundials_add_library(sundials_sunmatrixkokkosdense
    SOURCES
      sunmatrix_kokkosdense.cpp
    HEADERS
      ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_kokkosdense.hpp
    INCLUDE_SUBDIR
      sunmatrix
    OBJECT_LIBRARIES
      sundials_generic_obj
    LINK_LIBRARIES
      PUBLIC Kokkos::kokkos Kokkos::kokkoskernels
    OUTPUT_NAME
      sundials_sunmatrixkokkosdense
    VERSION
      ${sunmatrixlib_VERSION}
    SOVERSION
      ${sunmatrixlib_SOVERSION}
  )

  message(STATUS "Added SUNMATRIX_KOKKOSDENSE module")
endif()
