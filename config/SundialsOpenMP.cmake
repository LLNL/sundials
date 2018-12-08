# ---------------------------------------------------------------------------
# Programmer: David J. Gardner, and Cody J. Balos @ LLNL
# ---------------------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department
# of Energy by Lawrence Livermore National Laboratory in part under
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------------------
# Locate OpenMP and test for OpenMP device offloading support
# ---------------------------------------------------------------------------

find_package(OpenMP)

# Check for OpenMP offloading support
if(OPENMP_FOUND AND OPENMP_DEVICE_ENABLE)

  if(SKIP_OPENMP_DEVICE_CHECK)

    # The user has asked for checks to be skipped, assume offloading is supported
    set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)

  else()

    # If CMake version is 3.9 or newer, the FindOpenMP module checks the
    # OpenMP version.
    if((CMAKE_VERSION VERSION_EQUAL 3.9) OR (CMAKE_VERSION VERSION_GREATER 3.9))

      message(STATUS "Checking whether OpenMP supports device offloading")

      if((OpenMP_C_VERSION VERSION_EQUAL 4.5) OR (OpenMP_C_VERSION VERSION_GREATER 4.5))
        message(STATUS "Checking whether OpenMP supports device offloading -- yes")
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      else()
        message(STATUS "Checking whether OpenMP supports device offloading -- no")
        set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
      endif()

    else()

      # CMake OpenMP version check not available, assume offloading is supported
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
      print_warning("Unable to determine OpenMP offloading support."
        "OPENMP_DEVICE_ENABLE is ON but device offloading may not function.")

    endif()

  endif()

endif()
