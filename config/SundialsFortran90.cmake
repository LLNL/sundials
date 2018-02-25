# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2014, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# Fortran90-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

# If Fortran compiler flags are set using environemnt variables and both FFLAGS
# and FCFLAGS are used, then check if the variables are the same. If they are
# not the same then a fatal error occurs.
# 
# NOTE: This check must occur before 'enable_language(Fortran)' as it will use
# the value of FFLAGS to set CMAKE_Fortran_FLAGS
SET(ENV_FFLAGS "$ENV{FFLAGS}")
SET(ENV_FCFLAGS "$ENV{FCFLAGS}")
IF ((NOT "${ENV_FFLAGS}" STREQUAL "") AND
    (NOT "${ENV_FCFLAGS}" STREQUAL "") AND
    ("${CMAKE_Fortran_FLAGS}" STREQUAL ""))
  # check if environment variables are equal
  IF (NOT "${ENV_FFLAGS}" STREQUAL "${ENV_FCFLAGS}")
    PRINT_ERROR("FFLAGS='${ENV_FFLAGS}' and FCFLAGS='${ENV_FCFLAGS}' are both set but are not equal.")
  ENDIF()
ENDIF()

# Enable the language for next steps
enable_language(Fortran)
set(F90_FOUND TRUE)

# show some cache variables
MARK_AS_ADVANCED(CLEAR
  CMAKE_Fortran_COMPILER
  CMAKE_Fortran_FLAGS)

# hide all build type specific flags
MARK_AS_ADVANCED(FORCE
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_MINSIZEREL
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_RELWITHDEBINFO)

# only show flags for the current build type
# these flags are appended to CMAKE_Fortran_FLAGS
IF(CMAKE_BUILD_TYPE)       
  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE("Appending Fortran debug flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_DEBUG)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    MESSAGE("Appending Fortran min size release flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_MINSIZEREL)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE("Appending Fortran release flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_RELEASE)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    MESSAGE("Appending Fortran release with debug info flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
  ENDIF()
ENDIF()