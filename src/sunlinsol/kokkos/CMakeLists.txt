# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

if(BUILD_SUNLINSOL_KOKKOSDENSE)
  install(CODE "MESSAGE(\"\nInstall SUNLINSOL_KOKKOSDENSE\n\")")

  sundials_add_library(sundials_sunlinsolkokkosdense
    SOURCES
      sunlinsol_kokkosdense.cpp
    HEADERS
      ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_kokkosdense.hpp
    INCLUDE_SUBDIR
      sunlinsol
    OBJECT_LIBRARIES
      sundials_generic_obj
    LINK_LIBRARIES
      PUBLIC Kokkos::kokkos Kokkos::kokkoskernels
    OUTPUT_NAME
      sundials_sunlinsolkokkosdense
    VERSION
      ${sunlinsollib_VERSION}
    SOVERSION
      ${sunlinsollib_SOVERSION}
    )

  message(STATUS "Added SUNLINSOL_KOKKOSDENSE module")
endif()
