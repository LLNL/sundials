# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2017, Southern Methodist University and 
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
# CMakeLists.txt file for the band SUNMatrix library

INSTALL(CODE "MESSAGE(\"\nInstall SUNMATRIX_BAND\n\")")

# Add variable sunmatrixband_SOURCES with the sources for the SUNMATRIXBAND lib
SET(sunmatrixband_SOURCES sunmatrix_band.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXBAND library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunmatrixband_HEADERS with the exported SUNMATRIXBAND header files
SET(sunmatrixband_HEADERS sunmatrix_band.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunmatrix/ sunmatrixband_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXBAND library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunmatrixband_static STATIC ${sunmatrixband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixband_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixband CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunmatrixband_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXBAND library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunmatrixband_shared SHARED ${sunmatrixband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixband_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixband CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunmatrixband_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  INSTALL(TARGETS sundials_sunmatrixband_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNMATRIXBAND header files
INSTALL(FILES ${sunmatrixband_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXBAND library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunmatrixband_SOURCES fsunmatrix_band.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixband_static STATIC ${fsunmatrixband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixband_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixband CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunmatrixband_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixband_shared ${fsunmatrixband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixband CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunmatrixband_shared 
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    INSTALL(TARGETS sundials_fsunmatrixband_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNMATRIX_BAND module")
