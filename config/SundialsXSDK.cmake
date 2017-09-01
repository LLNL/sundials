# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# xSDK specific CMake variables. If set, these variables will 
# overwrite the corresponding SUNDIALS CMake variable.
#
# Only the USE_XSDK_DEFAULTS variable is created in CACHE by 
# default. The other variables are left undefined and can be set 
# by passing -D<variable_name>=<value> to cmake or they can be 
# activated in the cmake-gui by enabeling USE_XSDK_DEFAULTS.
# ---------------------------------------------------------------

# always show the option to turn on xSDK defaults
OPTION(USE_XSDK_DEFAULTS "Enable default xSDK settings" OFF)

IF(USE_XSDK_DEFAULTS)
  MESSAGE("Enabeling xSDK defaults")
  
  # set build type, SUNDIALS does not set a build type by default
  IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE("Setting build type to Debug")
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING 
      "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
  ENDIF()

  # set build precision, SUNDIALS_PRECISION defaults to double
  SET(XSDK_PRECISION "double" CACHE STRING "single, double, or quad")

  # set build index size, SUNDIALS_INDEX_TYPE defaults to int64_t
  SET(XSDK_INDEX_SIZE "32" CACHE STRING "32 or 64")

  # enable Fortran-C interface if at least one solver that provides such
  # an interface is built
  IF(BUILD_ARKODE OR BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)
    SHOW_VARIABLE(XSDK_ENABLE_FORTRAN BOOL "Enable Fortran-C support" OFF)
  ELSE()
    HIDE_VARIABLE(XSDK_ENABLE_FORTRAN BOOL "Enable Fortran-C support" OFF)
  ENDIF()

ENDIF()

