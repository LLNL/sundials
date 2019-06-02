# ---------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the MPIPlusX NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_MPIPLUSX\n\")")

IF(MPI_C_FOUND AND MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
ELSEIF(MPI_C_FOUND)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

# Add variable nvecmpiplusx_SOURCES with the sources for the NVECMPIPLUSX lib
SET(nvecmpiplusx_SOURCES
  nvector_mpiplusx.c
  ${sundials_SOURCE_DIR}/src/nvector/manyvector/nvector_manyvector.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECMPIPLUSX library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecmpiplusx_HEADERS with the exported NVECMPIPLUSX header files
SET(nvecmpiplusx_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_manyvector.h)
IF(MPI_ENABLE AND MPI_C_FOUND)
  LIST(APPEND nvecmpiplusx_HEADERS
      ${sundials_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h
      ${sundials_SOURCE_DIR}/include/nvector/nvector_mpiplusx.h)
ENDIF()

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECMPIPLUSX library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECMPIPLUSX library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecmpiplusx_static STATIC ${nvecmpiplusx_SOURCES} ${shared_SOURCES})

  SET_TARGET_PROPERTIES(sundials_nvecmpiplusx_static PROPERTIES
                        OUTPUT_NAME sundials_nvecmpiplusx
                        CLEAN_DIRECT_OUTPUT 1)

  INSTALL(TARGETS sundials_nvecmpiplusx_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECMPIPLUSX library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECMPIPLUSX library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecmpiplusx_shared SHARED ${nvecmpiplusx_SOURCES} ${shared_SOURCES})

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecmpiplusx_shared m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_nvecmpiplusx_shared PROPERTIES
                        OUTPUT_NAME sundials_nvecmpiplusx
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION}
                        CLEAN_DIRECT_OUTPUT 1)

  INSTALL(TARGETS sundials_nvecmpiplusx_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECMPIPLUSX header files
INSTALL(FILES ${nvecmpiplusx_HEADERS} DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_MPIPLUSX module")

