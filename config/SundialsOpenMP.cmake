# ---------------------------------------------------------------------------
# Programmer: David J. Gardner, and Cody J. Balos @ LLNL
# ---------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------------------
# Locate OpenMP and test for OpenMP verison and device offloading support
#
# Creates the variables:
#   OPENMP45_FOUND - was OpenMP v4.5 or greater found
#   OPENMP_SUPPORTS_DEVICE_OFFLOADING - is device offloading supported
# ---------------------------------------------------------------------------

set(OPENMP45_FOUND FALSE)
set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)

find_package(OpenMP)

# Work around a bug in setting OpenMP version variables in CMake >= 3.9. The
# OpenMP version information is not stored in cache variables and is not set
# on repeated calls to find OpenMP (i.e., when using ccmake). To ensure these
# variables exist store copies of the values.
set(OpenMP_C_VERSION "${OpenMP_C_VERSION}" CACHE INTERNAL "" FORCE)
set(OpenMP_CXX_VERSION "${OpenMP_CXX_VERSION}" CACHE INTERNAL "" FORCE)
set(OpenMP_Fortran_VERSION "${OpenMP_Fortran_VERSION}" CACHE INTERNAL "" FORCE)

# Check for OpenMP offloading support
if(OPENMP_FOUND AND (OPENMP_DEVICE_ENABLE OR SUPERLUDIST_OpenMP))

  if(SKIP_OPENMP_DEVICE_CHECK)

    # The user has asked for checks to be skipped, assume offloading is supported
    set(OPENMP45_FOUND TRUE)
    set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
    print_warning("Skipping OpenMP device/version check." "SUNDIALS OpenMP functionality dependent on OpenMP 4.5+ is not guaranteed.")

  else()

    # If CMake version is 3.9 or newer, the FindOpenMP module checks the OpenMP version.
    if((CMAKE_VERSION VERSION_EQUAL 3.9) OR (CMAKE_VERSION VERSION_GREATER 3.9))

      message(STATUS "Checking whether OpenMP supports device offloading")
      if((OpenMP_C_VERSION VERSION_EQUAL 4.5) OR (OpenMP_C_VERSION VERSION_GREATER 4.5))
        message(STATUS "Checking whether OpenMP supports device offloading -- yes")
        set(OPENMP45_FOUND TRUE)
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      else()
        message(STATUS "Checking whether OpenMP supports device offloading -- no")
        set(OPENMP45_FOUND FALSE)
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
      endif()

    else()

      # CMake OpenMP version check not available. Assume 4.5+ and that offloading is supported.
      set(OPENMP45_FOUND TRUE)
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      print_warning("Unable to determine OpenMP offloading support." "SUNDIALS OpenMP functionality dependent on OpenMP 4.5+ is not guaranteed.")

    endif()

  endif()

endif()
