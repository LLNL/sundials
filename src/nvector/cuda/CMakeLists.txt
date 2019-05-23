# ---------------------------------------------------------------
# Programmer:  Slaven Peles, and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_CUDA\n\")")

# Add variable nveccuda_SOURCES with the sources for the NVECUDA lib
SET(nveccuda_SOURCES nvector_cuda.cu)

# Tell compiler it is a CUDA source
set_source_files_properties(${nveccuda_SOURCES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECCUDA library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECUDA header files
SET(nveccuda_HEADERS nvector_cuda.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nveccuda_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECCUDA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECUDA library
IF(BUILD_STATIC_LIBS)
  CUDA_ADD_LIBRARY(sundials_nveccuda_static STATIC ${nveccuda_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nveccuda_static PROPERTIES
                        OUTPUT_NAME sundials_nveccuda
                        CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nveccuda_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECUDA library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECUDA library
IF(BUILD_SHARED_LIBS)
  CUDA_ADD_LIBRARY(sundials_nveccuda_shared SHARED ${nveccuda_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nveccuda_shared PROPERTIES
                        OUTPUT_NAME sundials_nveccuda
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nveccuda_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the CUDA NVector header files
INSTALL(FILES ${nveccuda_HEADERS} DESTINATION include/nvector)
INSTALL(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/cuda DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_CUDA module")
