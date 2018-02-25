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
# C++-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

# Enable the language for next steps
enable_language(CXX)
set(CXX_FOUND TRUE)

# show some cache variables
MARK_AS_ADVANCED(CLEAR
  CMAKE_CXX_COMPILER
  CMAKE_CXX_FLAGS)

# hide all build type specific flags
MARK_AS_ADVANCED(FORCE
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_CXX_FLAGS_RELWITHDEBINFO)

# only show flags for the current build type
# these flags are appended to CMAKE_CXX_FLAGS
IF(CMAKE_BUILD_TYPE)
  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE("Appending CXX debug flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_DEBUG)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    MESSAGE("Appending CXX min size release flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_MINSIZEREL)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE("Appending CXX release flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELEASE)
  ELSEIF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    MESSAGE("Appending CXX release with debug info flags")
    MARK_AS_ADVANCED(CLEAR CMAKE_CXX_FLAGS_RELWITHDEBINFO)
  ENDIF()
ENDIF()