# ---------------------------------------------------------------
# Programmer:  Carol S. Woodward
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the spc_parallel NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_SPCPAR\n\")")

IF(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSE()
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable nvecspcpar_SOURCES with the sources for the NVECSPCPAR lib
SET(nvecspcpar_SOURCES nvector_spcparallel.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSPCPAR library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecspcpar_HEADERS with the exported NVECSPCPAR header files
SET(nvecspcpar_HEADERS nvector_spcparallel.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecspcpar_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECSPCPAR library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSPCPAR library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecspcpar_static STATIC ${nvecspcpar_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecspcpar_static
    PROPERTIES OUTPUT_NAME sundials_nvecspcpar CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecspcpar_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSPCPAR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSPCPAR library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecspcpar_shared SHARED ${nvecspcpar_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecspcpar_shared
    PROPERTIES OUTPUT_NAME sundials_nvecspcpar CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecspcpar_shared 
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecspcpar_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECSPCPAR header files
INSTALL(FILES ${nvecspcpar_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI works, build and install the FNVECSPCPAR library
IF(F77_INTERFACE_ENABLE AND MPI_Fortran_FOUND)
  SET(fnvecspcpar_SOURCES fnvector_spcparallel.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecspcpar_static STATIC ${fnvecspcpar_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecspcpar_static
      PROPERTIES OUTPUT_NAME sundials_fnvecspcpar CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecspcpar_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecspcpar_shared SHARED ${fnvecspcpar_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecspcpar_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecspcpar CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecspcpar_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecspcpar_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF(BUILD_SHARED_LIBS)
ENDIF()
#
MESSAGE(STATUS "Added NVECTOR_SPCPAR module")
