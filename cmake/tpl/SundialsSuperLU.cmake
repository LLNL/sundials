# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup SUPERLU correctly.
# Created from the SundialsTPL.cmake.template template.
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

if(NOT DEFINED SUNDIALS_SUPERLU_INCLUDED)
  set(SUNDIALS_SUPERLU_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# SUPERLU does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  print_error("SUPERLU is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

if(NOT (SUNDIALS_INDEX_SIZE MATCHES "32"))
  print_error("SuperLU requires SUNDIALS_INDEX_SIZE=32")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(SUPERLU REQUIRED)

message(STATUS "SUPERLU_LIBRARIES:   ${SUPERLU_LIBRARIES}")
message(STATUS "SUPERLU_INCLUDE_DIR: ${SUPERLU_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUPERLU_FOUND AND (NOT SUPERLU_WORKS))

  # Create the SUPERLU_TEST directory
  set(SUPERLU_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLU_TEST)
  file(MAKE_DIRECTORY ${SUPERLU_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${SUPERLU_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.1.3)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_INCLUDE_DIRECTORIES(ltest PRIVATE ${SUPERLU_INCLUDE_DIR})\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLU_LIBRARIES})\n")

  # Create a C source file which calls a SUPERLU function
  file(WRITE ${SUPERLU_TEST_DIR}/ltest.c
    "\#include \"slu_ddefs.h\"\n"
    "int main(){\n"
    "SuperMatrix *A;\n"
    "NCformat *Astore;\n"
    "A = NULL;\n"
    "Astore = NULL;\n"
    "if (A != NULL || Astore != NULL) return(1);\n"
    "else return(0);\n"
    "}\n")


  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${SUPERLU_TEST_DIR} ${SUPERLU_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${SUPERLU_TEST_DIR}/CMakeFiles)

 # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if SuperLU works with SUNDIALS... OK")
    set(SUPERLU_WORKS TRUE CACHE BOOL "SuperLU works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if SuperLU works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to SuperLU is not functional.")
  endif()

elseif(SUPERLU_FOUND AND SUPERLU_WORKS)
  message(STATUS "Skipped SuperLU tests, assuming SuperLU works with SUNDIALS. Set SUPERLU_WORKS=FALSE to (re)run compatibility test.")
endif()
