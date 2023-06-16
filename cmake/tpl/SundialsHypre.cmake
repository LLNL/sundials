# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2023, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup HYPRE correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_HYPRE_INCLUDED)
  set(SUNDIALS_HYPRE_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

if(ENABLE_HYPRE)
  # Using hypre requires building with MPI enabled
  if(NOT ENABLE_MPI)
    print_error("MPI is required for hypre support. Set ENABLE_MPI to ON.")
  endif()
  # Using hypre requres C99 or newer
  if(CMAKE_C_STANDARD STREQUAL "90")
    message(SEND_ERROR "CMAKE_C_STANDARD must be >= c99 with ENABLE_HYPRE=ON")
  endif()
endif()

if((SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA") AND (NOT ENABLE_CUDA))
  message(FATAL_ERROR "hypre with a CUDA backend requires ENABLE_CUDA = ON")
endif()

if((SUNDIALS_HYPRE_BACKENDS MATCHES "HIP") AND (NOT ENABLE_HIP))
  message(FATAL_ERROR "hypre with a HIP backend requires ENABLE_HIP = ON")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(HYPRE REQUIRED)

message(STATUS "HYPRE_LIBRARIES:   ${HYPRE_LIBRARIES}")
message(STATUS "HYPRE_INCLUDE_DIR: ${HYPRE_INCLUDE_DIR}")

#TODO verify this is good code -jsdomine

# Find the hypre library configuration file (HYPRE_config.h)
find_file(HYPRE_CONFIGH_PATH HYPRE_config.h
          HINTS "${HYPRE_DIR}"
          PATH_SUFFIXES include
          NO_DEFAULT_PATH
          REQUIRED) #TODO should we set as required? -jsdomine
mark_as_advanced(FORCE HYPRE_CONFIGH_PATH)

# --- Throws error if hypre wasn't built with CMake; HYPREConfig.cmake is not used, anyway ---
# # Look for CMake configuration file in hypre installation (HYPREConfig.cmake)
# find_package(HYPRE CONFIG
#              HINTS "${HYPRE_DIR}"
#              NO_DEFAULT_PATH
#              REQUIRED)

# Display hypre version
file(READ "${HYPRE_CONFIGH_PATH}" _hypre_config_file_text)
string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" _hypre_release_version "${_hypre_config_file_text}")
message(STATUS "hypre Version: ${_hypre_release_version}")

# Determine the backends
foreach(_backend CUDA HIP)
  file(STRINGS "${HYPRE_CONFIGH_PATH}" _hypre_has_backend REGEX "^#define HYPRE_USING_${_backend}")
  if(_hypre_has_backend)
    set(HYPRE_BACKENDS "${_backend};${HYPRE_BACKENDS}")
    message(STATUS "hypre built with ${_backend} backend? - YES")
  else()
    message(STATUS "hypre built with ${_backend} backend? - NO")
  endif()
endforeach()

# Check for CUDA Unified Memory
file(STRINGS "${HYPRE_CONFIGH_PATH}" _hypre_using_unified_memory REGEX "^#define HYPRE_USING_UNIFIED_MEMORY")
if(_hypre_using_unified_memory)
  set(SUNDIALS_HYPRE_USING_UNIFIED_MEMORY TRUE)
  message(STATUS "hypre using CUDA Unified Memory? - YES")
else()
  message(STATUS "hypre using CUDA Unified Memory? - NO")
endif()

# Manually link to cuda_runtime (CUDA::cudart) when using CUDA backend
# set(HYPRE_NEEDS_THREADS OFF)
if(SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
  list(APPEND HYPRE_LIBRARIES CUDA::cudart CUDA::cublas)
#   set(HYPRE_NEEDS_THREADS ON)
#   if(NOT TARGET Threads::Threads)
#     find_package(Threads)
#   endif()
#   if(NOT TARGET cuda_runtime)
#     add_library(cuda_runtime INTERFACE IMPORTED)
#     target_link_libraries(cuda_runtime INTERFACE CUDA::cudart)
#   endif()
endif()

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

message(STATUS "Requested SUNDIALS hypre backend: ${SUNDIALS_HYPRE_BACKENDS}")

if((SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA") AND
   (NOT HYPRE_BACKENDS MATCHES "CUDA"))
  print_error("Requested that SUNDIALS uses the CUDA hypre backend, but hypre was not built with the CUDA backend.")
endif()

if((SUNDIALS_HYPRE_BACKENDS MATCHES "HIP") AND
   (NOT HYPRE_BACKENDS MATCHES "HIP"))
  print_error("Requested that SUNDIALS uses the HIP hypre backend, but hypre was not built with the HIP backend.")
endif()

if(HYPRE_FOUND AND (NOT HYPRE_WORKS))
  # Do any checks which don't require compilation first.

  # Create the HYPRE_TEST directory
  set(HYPRE_TEST_DIR ${PROJECT_BINARY_DIR}/HYPRE_TEST)
  file(MAKE_DIRECTORY ${HYPRE_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${HYPRE_TEST_DIR}/CMakeLists.txt
  "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
  "PROJECT(ltest C)\n"
  "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
  "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
  "SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})\n"
  "SET(CMAKE_C_STANDARD \"${CMAKE_C_STANDARD}\")\n"
  "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
  "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
  "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
  "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
  "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
  "SET(CMAKE_EXE_LINKER_FLAGS \"${LINK_MATH_LIB}\")\n"
  "INCLUDE_DIRECTORIES(${HYPRE_INCLUDE_DIR})\n"
  "ADD_EXECUTABLE(ltest ltest.c)\n"
  "TARGET_LINK_LIBRARIES(ltest ${HYPRE_LIBRARIES})\n")

  file(WRITE ${HYPRE_TEST_DIR}/ltest.c
  "\#include \"HYPRE_parcsr_ls.h\"\n"
  "int main(){\n"
  "HYPRE_ParVector par_b;\n"
  "HYPRE_IJVector b;\n"
  "par_b = 0;\n"
  "b = 0;\n"
  "if (par_b != 0 || b != 0) return(1);\n"
  "else return(0);\n"
  "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${HYPRE_TEST_DIR}/CMakeFiles)

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${HYPRE_TEST_DIR} ${HYPRE_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if hypre works... OK")
    set(HYPRE_WORKS TRUE CACHE BOOL "hypre works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if hypre works... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to hypre is not functional.")
  endif()

elseif(HYPRE_FOUND AND HYPRE_WORKS)
  message(STATUS "Skipped hypre tests, assuming hypre works with SUNDIALS. Set HYPRE_WORKS=FALSE to (re)run compatibility test.")
endif()
