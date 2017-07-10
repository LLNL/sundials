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
# CMakeLists.txt file for the SPFGMR SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPFGMR\n\")")

# Add variable sunlinsolspfgmr_SOURCES with the sources for the SUNLINSOLSPFGMR lib
SET(sunlinsolspfgmr_SOURCES sunlinsol_spfgmr.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPFGMR library
SET(shared_SOURCES
  sundials_math.c
  sundials_iterative.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunlinsolspfgmr_HEADERS with the exported SUNLINSOLSPFGMR header files
SET(sunlinsolspfgmr_HEADERS sunlinsol_spfgmr.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunlinsol/ sunlinsolspfgmr_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPFGMR library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPFGMR library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspfgmr_static STATIC ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsolspfgmr_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLSPFGMR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLSPFGMR library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsolspfgmr_shared SHARED ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsolspfgmr_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLSPFGMR header files
INSTALL(FILES ${sunlinsolspfgmr_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLSPFGMR library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsolspfgmr_SOURCES fsunlinsol_spfgmr.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspfgmr_static STATIC ${fsunlinsolspfgmr_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsolspfgmr_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsolspfgmr_shared ${fsunlinsolspfgmr_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsolspfgmr_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsolspfgmr_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_SPFGMR module")
