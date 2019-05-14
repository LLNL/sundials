# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the ManyVector NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_MANYVECTOR\n\")")

# Add variable nvecmanyvector_SOURCES with the sources for the NVECMANYVECTOR lib
SET(nvecmanyvector_SOURCES nvector_manyvector.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECMANYVECTOR library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecmanyvector_HEADERS with the exported NVECMANYVECTOR header files
SET(nvecmanyvector_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_manyvector.h)

IF(MPI_ENABLE AND MPI_C_FOUND)
  INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_MPIMANYVECTOR\n\")")
  LIST(APPEND nvecmanyvector_HEADERS ${sundials_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h)
  IF(MPI_C_COMPILER)
    # use MPI wrapper as the compiler
    SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  ELSEIF()
    # add MPI_INCLUDE_PATH to include directories
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  ENDIF()
ENDIF()

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECMANYVECTOR library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECMANYVECTOR library
IF(BUILD_STATIC_LIBS)

  # ----------------------- No MPI
  ADD_LIBRARY(sundials_nvecmanyvector_static STATIC ${nvecmanyvector_SOURCES} ${shared_SOURCES})

  TARGET_COMPILE_DEFINITIONS(sundials_nvecmanyvector_static PRIVATE -DSUNDIALS_MPI_ENABLED=0)

  SET_TARGET_PROPERTIES(sundials_nvecmanyvector_static PROPERTIES
                        OUTPUT_NAME sundials_nvecmanyvector
                        CLEAN_DIRECT_OUTPUT 1)

  INSTALL(TARGETS sundials_nvecmanyvector_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_LIBRARY(sundials_nvecmpimanyvector_static STATIC ${nvecmanyvector_SOURCES} ${shared_SOURCES})

    TARGET_COMPILE_DEFINITIONS(sundials_nvecmpimanyvector_static PRIVATE -DSUNDIALS_MPI_ENABLED=1)

    SET_TARGET_PROPERTIES(sundials_nvecmpimanyvector_static PROPERTIES
                          OUTPUT_NAME sundials_nvecmpimanyvector
                          CLEAN_DIRECT_OUTPUT 1)

    INSTALL(TARGETS sundials_nvecmpimanyvector_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()

ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECMANYVECTOR library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECMANYVECTOR library
IF(BUILD_SHARED_LIBS)

  # ----------------------- No MPI
  ADD_LIBRARY(sundials_nvecmanyvector_shared SHARED ${nvecmanyvector_SOURCES} ${shared_SOURCES})

  TARGET_COMPILE_DEFINITIONS(sundials_nvecmanyvector_shared PRIVATE -DSUNDIALS_MPI_ENABLED=0)

  IF(UNIX)
    TARGET_LINK_LIBRARIES(sundials_nvecmanyvector_shared PRIVATE m)
  ENDIF()

  SET_TARGET_PROPERTIES(sundials_nvecmanyvector_shared PROPERTIES
                        OUTPUT_NAME sundials_nvecmanyvector
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION}
                        CLEAN_DIRECT_OUTPUT 1)

  INSTALL(TARGETS sundials_nvecmanyvector_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})

  # ----------------------- With MPI
  IF(MPI_ENABLE AND MPI_C_FOUND)
    ADD_LIBRARY(sundials_nvecmpimanyvector_shared SHARED ${nvecmanyvector_SOURCES} ${shared_SOURCES})

    TARGET_COMPILE_DEFINITIONS(sundials_nvecmpimanyvector_shared PRIVATE -DSUNDIALS_MPI_ENABLED=1)

    IF(UNIX)
      TARGET_LINK_LIBRARIES(sundials_nvecmpimanyvector_shared PRIVATE m)
    ENDIF()

    SET_TARGET_PROPERTIES(sundials_nvecmpimanyvector_shared PROPERTIES
                          OUTPUT_NAME sundials_nvecmpimanyvector
                          VERSION ${nveclib_VERSION}
                          SOVERSION ${nveclib_SOVERSION}
                          CLEAN_DIRECT_OUTPUT 1)

    INSTALL(TARGETS sundials_nvecmpimanyvector_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
  ENDIF()

ENDIF(BUILD_SHARED_LIBS)

# Install the NVECMANYVECTOR header files
INSTALL(FILES ${nvecmanyvector_HEADERS} DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_MANYVECTOR module")

IF(MPI_ENABLE AND MPI_C_FOUND)
  MESSAGE(STATUS "Added NVECTOR_MPIMANYVECTOR module")
ENDIF()

