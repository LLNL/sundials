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

# --- hypre requires MPI and C99 or newer ---
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

# --- hypre with CUDA backend requires bulding with CUDA enabled ---
if((SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA") AND (NOT ENABLE_CUDA))
  message(FATAL_ERROR "hypre with a CUDA backend requires ENABLE_CUDA = ON")
endif()

# --- hypre with HIP backend requires bulding with HIP enabled ---
if((SUNDIALS_HYPRE_BACKENDS MATCHES "HIP") AND (NOT ENABLE_HIP))
  message(FATAL_ERROR "hypre with a HIP backend requires ENABLE_HIP = ON")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(HYPRE REQUIRED)

message(STATUS "HYPRE_LIBRARIES:   ${HYPRE_LIBRARIES}")
message(STATUS "HYPRE_INCLUDE_DIR: ${HYPRE_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# --- Ensure hypre built with requested backend ---
message(STATUS "Requested SUNDIALS hypre backend: ${SUNDIALS_HYPRE_BACKENDS}")
if((SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA") AND
   (NOT HYPRE_BACKENDS MATCHES "CUDA"))
  print_error("Requested that SUNDIALS uses the CUDA hypre backend, but hypre was not built with the CUDA backend.")
endif()
if((SUNDIALS_HYPRE_BACKENDS MATCHES "HIP") AND
   (NOT HYPRE_BACKENDS MATCHES "HIP"))
  print_error("Requested that SUNDIALS uses the HIP hypre backend, but hypre was not built with the HIP backend.")
endif()

# --- Test hypre libraries ---
if(HYPRE_FOUND AND (NOT HYPRE_WORKS))
  # Do any checks which don't require compilation first.

  # Create the HYPRE_TEST directory
  set(HYPRE_TEST_DIR ${PROJECT_BINARY_DIR}/HYPRE_TEST)
  file(MAKE_DIRECTORY ${HYPRE_TEST_DIR})

  # Choose test language (C for serial, C++ for CUDA or HIP)
  if (SUNDIALS_HYPRE_BACKENDS MATCHES "SERIAL")
    set(_hypre_test_project_languages C)
    set(_hypre_test_file_language C)
  elseif (SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
    set(_hypre_test_project_languages C CXX CUDA)
    set(_hypre_test_file_language CUDA)
  elseif (SUNDIALS_HYPRE_BACKENDS MATCHES "HIP")
    set(_hypre_test_project_languages C CXX HIP)
    set(_hypre_test_file_language HIP)
  endif ()

  # Create a CMakeLists.txt file
  file(WRITE ${HYPRE_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ltest ${_hypre_test_project_languages})\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n")
  foreach(_language IN LISTS _hypre_test_project_languages)
    file(APPEND ${HYPRE_TEST_DIR}/CMakeLists.txt
      "SET(CMAKE_${_language}_COMPILER ${MPI_${_language}_COMPILER})\n"
      "SET(CMAKE_${_language}_STANDARD \"${CMAKE_${_language}_STANDARD}\")\n"
      "SET(CMAKE_${_language}_FLAGS \"${CMAKE_${_language}_FLAGS}\")\n"
      "SET(CMAKE_${_language}_FLAGS_RELEASE \"${CMAKE_${_language}_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_${_language}_FLAGS_DEBUG \"${CMAKE_${_language}_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_${_language}_FLAGS_RELWITHDEBUGINFO \"${CMAKE_${_language}_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_${_language}_FLAGS_MINSIZE \"${CMAKE_${_language}_FLAGS_MINSIZE}\")\n")
  endforeach()
  file(APPEND ${HYPRE_TEST_DIR}/CMakeLists.txt
    "SET(CMAKE_EXE_LINKER_FLAGS \"${LINK_MATH_LIB}\")\n"
    "INCLUDE_DIRECTORIES(${HYPRE_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${HYPRE_LIBRARIES})\n"
    "SET_SOURCE_FILES_PROPERTIES(ltest.c LANGUAGE ${_hypre_test_file_language})\n")

  if(SUNDIALS_HYPRE_BACKENDS MATCHES "CUDA")
    file(APPEND ${HYPRE_TEST_DIR}/CMakeLists.txt
      "TARGET_COMPILE_FEATURES(ltest PUBLIC USE_CUDA)\n")
  endif()

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
