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
# CMakeLists.txt file for the diagonal SUNMatrix library

INSTALL(CODE "MESSAGE(\"\nInstall SUNMATRIX_DIAGONAL\n\")")

# Add variable sunmatrixdiagonal_SOURCES with the sources for the SUNMATRIXDIAGONAL lib
SET(sunmatrixdiagonal_SOURCES sunmatrix_diagonal.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNMATRIXDIAGONAL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunmatrixdiagonal_HEADERS with the exported SUNMATRIXDIAGONAL header files
SET(sunmatrixdiagonal_HEADERS sunmatrix_diagonal.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunmatrix/ sunmatrixdiagonal_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNMATRIXDIAGONAL library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNMATRIXDIAGONAL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunmatrixdiagonal_static STATIC ${sunmatrixdiagonal_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixdiagonal_static
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdiagonal CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunmatrixdiagonal_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNMATRIXDIAGONAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNMATRIXDIAGONAL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunmatrixdiagonal_shared SHARED ${sunmatrixdiagonal_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunmatrixdiagonal_shared
    PROPERTIES OUTPUT_NAME sundials_sunmatrixdiagonal CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunmatrixdiagonal_shared
    PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
  INSTALL(TARGETS sundials_sunmatrixdiagonal_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNMATRIXDIAGONAL header files
INSTALL(FILES ${sunmatrixdiagonal_HEADERS} DESTINATION include/sunmatrix)

# If FCMIX is enabled, build and install the FSUNMATRIXDIAGONAL library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunmatrixdiagonal_SOURCES fsunmatrix_diagonal.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixdiagonal_static STATIC ${fsunmatrixdiagonal_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixdiagonal_static
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdiagonal CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunmatrixdiagonal_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunmatrixdiagonal_shared ${fsunmatrixdiagonal_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunmatrixdiagonal_shared
      PROPERTIES OUTPUT_NAME sundials_fsunmatrixdiagonal CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunmatrixdiagonal_shared 
      PROPERTIES VERSION ${sunmatrixlib_VERSION} SOVERSION ${sunmatrixlib_SOVERSION})
    INSTALL(TARGETS sundials_fsunmatrixdiagonal_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNMATRIX_DIAGONAL module")
