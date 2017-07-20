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
# CMakeLists.txt file for the LapackBand SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKBAND\n\")")

# Add variable sunlinsollapackband_SOURCES with the sources for the SUNLINSOLLAPACKBAND lib
SET(sunlinsollapackband_SOURCES sunlinsol_lapackband.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLLAPACKBAND library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunlinsollapackband_HEADERS with the exported SUNLINSOLLAPACKBAND header files
SET(sunlinsollapackband_HEADERS sunlinsol_lapackband.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunlinsol/ sunlinsollapackband_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLLAPACKBAND library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLLAPACKBAND library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackband_static STATIC ${sunlinsollapackband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsollapackband_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLLAPACKBAND library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLLAPACKBAND library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsollapackband_shared SHARED ${sunlinsollapackband_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsollapackband_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsollapackband_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLLAPACKBAND header files
INSTALL(FILES ${sunlinsollapackband_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLLAPACKBAND library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsollapackband_SOURCES fsunlinsol_lapackband.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackband_static STATIC ${fsunlinsollapackband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsollapackband_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsollapackband_shared ${fsunlinsollapackband_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsollapackband CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsollapackband_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsollapackband_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_LAPACKBAND module")
