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
# Module to find and setup CALIPER correctly.
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

if(NOT DEFINED SUNDIALS_CALIPER_INCLUDED)
  set(SUNDIALS_CALIPER_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(CALIPER
             PATHS "${CALIPER_DIR}"
             REQUIRED)

message(STATUS "CALIPER_LIB_DIR:     ${caliper_LIB_DIR}")
message(STATUS "CALIPER_INCLUDE_DIR: ${caliper_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(CALIPER_FOUND AND (NOT CALIPER_WORKS))
  # Do any checks which don't require compilation first.

  # Create the CALIPER_TEST directory
  set(CALIPER_TEST_DIR ${PROJECT_BINARY_DIR}/CALIPER_TEST)
  file(MAKE_DIRECTORY ${CALIPER_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${CALIPER_TEST_DIR}/CMakeLists.txt
    "cmake_minimum_required(VERSION ${CMAKE_VERSION})\n"
    "project(ltest C)\n"
    "set(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "set(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "set(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "set(CMAKE_C_STANDARD \"${CMAKE_C_STANDARD}\")\n"
    "set(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "set(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "set(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "set(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "set(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "add_executable(ltest ltest.c)\n"
    "target_include_directories(ltest PRIVATE \"${caliper_INCLUDE_DIR}\")\n"
    "target_link_libraries(ltest \"-L${caliper_LIB_DIR}\")\n"
    "target_link_libraries(ltest caliper)\n")

  # Create a C source file
  file(WRITE ${CALIPER_TEST_DIR}/ltest.c
  "\#include <caliper/cali.h>\n"
  "int main()\n"
  "{\n"
  "  CALI_MARK_FUNCTION_BEGIN;\n"
  "  CALI_MARK_FUNCTION_END;\n"
  "  return 0;\n"
  "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${CALIPER_TEST_DIR}/CMakeFiles)

  # Attempt to build and link the "ltest" executable
  try_compile(COMPILE_OK ${CALIPER_TEST_DIR} ${CALIPER_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if CALIPER works with SUNDIALS... OK")
    set(CALIPER_WORKS TRUE CACHE BOOL "CALIPER works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if CALIPER works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    print_error("SUNDIALS interface to CALIPER is not functional.")
  endif()

elseif(CALIPER_FOUND AND CALIPER_WORKS)
  message(STATUS "Skipped CALIPER tests, assuming CALIPER works with SUNDIALS.")
endif()
