# ---------------------------------------------------------------
# Programmer:  Slaven Peles, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the RAJA NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_RAJA\n\")")

# Add variable nvecraja_SOURCES with the sources for the NVECRAJA lib
SET(nvecraja_SOURCES
  nvector_raja.cu
)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECRAJA library
SET(shared_SOURCES
  sundials_nvector.c
  sundials_math.c
)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nveccuda_HEADERS with the exported NVECRAJA header files
SET(nvecraja_HEADERS nvector_raja.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecraja_HEADERS)

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECRAJA library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nveccudaraja_static STATIC ${nvecraja_SOURCES} ${shared_SOURCES})
  TARGET_COMPILE_FEATURES(sundials_nveccudaraja_static PUBLIC cxx_std_11)
  TARGET_COMPILE_OPTIONS(sundials_nveccudaraja_static PUBLIC
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>)
  TARGET_LINK_LIBRARIES(sundials_nveccudaraja_static PUBLIC RAJA)
  SET_TARGET_PROPERTIES(sundials_nveccudaraja_static PROPERTIES
                        OUTPUT_NAME sundials_nveccudaraja
                        CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nveccudaraja_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECRAJA library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECRAJA library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nveccudaraja_shared SHARED ${nvecraja_SOURCES} ${shared_SOURCES})
  TARGET_COMPILE_FEATURES(sundials_nveccudaraja_shared PUBLIC cxx_std_11)
  TARGET_COMPILE_OPTIONS(sundials_nveccudaraja_shared PUBLIC
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>)
  TARGET_LINK_LIBRARIES(sundials_nveccudaraja_shared PUBLIC RAJA)
  SET_TARGET_PROPERTIES(sundials_nveccudaraja_shared PROPERTIES
                        OUTPUT_NAME sundials_nveccudaraja
                        CLEAN_DIRECT_OUTPUT 1
                        VERSION ${nveclib_VERSION}
                        SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nveccudaraja_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECRAJA header files
INSTALL(FILES ${nvecraja_HEADERS} DESTINATION include/nvector)
INSTALL(DIRECTORY ${sundials_SOURCE_DIR}/include/nvector/raja DESTINATION include/nvector)

MESSAGE(STATUS "Added NVECTOR_RAJA module")
