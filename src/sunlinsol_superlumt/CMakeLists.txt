# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds, Ashley Crawford @ SMU
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
# CMakeLists.txt file for the band SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_BAND\n\")")

# Add variable sunlinsolband_SOURCES with the sources for the SUNLINSOLBAND lib
SET(sunlinsolband_SOURCES sunlinsol_band.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLBAND library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunlinsolband_HEADERS with the exported SUNLINSOLBAND header files
SET(sunlinsolband_HEADERS sunlinsol_band.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunlinsol/ sunlinsolband_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLBAND library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolband_static STATIC ${sunlinsolband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolband_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolband CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolband_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLBAND library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolband_shared SHARED ${sunlinsolband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolband_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolband CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolband_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolband_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLBAND header files
INSTALL(FILES ${sunlinsolband_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLBAND library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolband_SOURCES fsunlinsol_band.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolband_static STATIC ${fsunlinsolband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolband CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolband_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolband_shared ${fsunlinsolband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolband CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolband_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolband_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_BAND module")
