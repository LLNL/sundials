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
# CMakeLists.txt file for the diagonal SUNLinearSolver library

INSTALL(CODE "MESSAGE(\"\nInstall SUNLINSOL_DIAGONAL\n\")")

# Add variable sunlinsoldiagonal_SOURCES with the sources for the SUNLINSOLDIAGONAL lib
SET(sunlinsoldiagonal_SOURCES sunlinsol_diagonal.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLDIAGONAL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable sunlinsoldiagonal_HEADERS with the exported SUNLINSOLDIAGONAL header files
SET(sunlinsoldiagonal_HEADERS sunlinsol_diagonal.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/sunlinsol/ sunlinsoldiagonal_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLDIAGONAL library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLDIAGONAL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_sunlinsoldiagonal_static STATIC ${sunlinsoldiagonal_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsoldiagonal_static
    PROPERTIES OUTPUT_NAME sundials_sunlinsoldiagonal CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_sunlinsoldiagonal_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the SUNLINSOLDIAGONAL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the SUNLINSOLDIAGONAL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_sunlinsoldiagonal_shared SHARED ${sunlinsoldiagonal_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_sunlinsoldiagonal_shared
    PROPERTIES OUTPUT_NAME sundials_sunlinsoldiagonal CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_sunlinsoldiagonal_shared
    PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
  INSTALL(TARGETS sundials_sunlinsoldiagonal_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the SUNLINSOLDIAGONAL header files
INSTALL(FILES ${sunlinsoldiagonal_HEADERS} DESTINATION include/sunlinsol)

# If FCMIX is enabled, build and install the FSUNLINSOLDIAGONAL library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fsunlinsoldiagonal_SOURCES fsunlinsol_diagonal.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fsunlinsoldiagonal_static STATIC ${fsunlinsoldiagonal_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsoldiagonal_static
      PROPERTIES OUTPUT_NAME sundials_fsunlinsoldiagonal CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fsunlinsoldiagonal_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fsunlinsoldiagonal_shared ${fsunlinsoldiagonal_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fsunlinsoldiagonal_shared
      PROPERTIES OUTPUT_NAME sundials_fsunlinsoldiagonal CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fsunlinsoldiagonal_shared 
      PROPERTIES VERSION ${sunlinsollib_VERSION} SOVERSION ${sunlinsollib_SOVERSION})
    INSTALL(TARGETS sundials_fsunlinsoldiagonal_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added SUNLINSOL_DIAGONAL module")
