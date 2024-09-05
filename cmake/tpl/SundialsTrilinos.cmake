# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup Trilinos correctly.
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

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# Find Trilinos
find_package(Trilinos REQUIRED COMPONENTS Tpetra
  HINTS "${Trilinos_DIR}/lib/cmake/Trilinos" "${Trilinos_DIR}")

message(STATUS "Trilinos_LIBRARIES:    ${Trilinos_LIBRARIES}")
message(STATUS "Trilinos_INCLUDE_DIRS: ${Trilinos_INCLUDE_DIRS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Try building a simple test
if(NOT Trilinos_WORKS)

  message(CHECK_START "Testing Trilinos")

  # Create the test directory
  set(TRILINOS_TEST_DIR ${PROJECT_BINARY_DIR}/TRILINOS_TEST)

  # Create a CXX source file calling a Trilinos Tpetra function
  file(
    WRITE ${TRILINOS_TEST_DIR}/test.cxx
    "#include <Tpetra_Version.hpp>\n"
    "int main(void) {\n"
    "std::cout << Tpetra::version() << std::endl;\n"
    "return(0);\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TRILINOS_TEST_DIR}
    ${TRILINOS_TEST_DIR}/test.cxx
    LINK_LIBRARIES Tpetra::all_libs
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check the result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TRILINOS_TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
      "Could not compile Trilinos test. Check output in ${TRILINOS_TEST_DIR}/compile.out"
    )
  endif()

else()
  message(STATUS "Skipped Trilinos test. Set TRILINOS_WORKS=FALSE to test.")
endif()
