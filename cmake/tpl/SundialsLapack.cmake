# -----------------------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module to find and setup LAPACK/BLAS correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# <TPL>_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# LAPACK does not support extended precision
if(ENABLE_LAPACK AND SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "LAPACK is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(LAPACK REQUIRED)

# get path to LAPACK library to use in generated makefiles for examples, if
# LAPACK_LIBRARIES contains multiple items only use the path of the first entry
list(GET LAPACK_LIBRARIES 0 TMP_LAPACK_LIBRARIES)
get_filename_component(LAPACK_LIBRARY_DIR ${TMP_LAPACK_LIBRARIES} PATH)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Macro to test different name mangling options
macro(test_lapack_name_mangling)

  # Both the name mangling case and number of underscores must be set
  if((NOT SUNDIALS_LAPACK_FUNC_CASE AND SUNDIALS_LAPACK_FUNC_UNDERSCORES) OR
      (SUNDIALS_LAPACK_FUNC_CASE AND NOT SUNDIALS_LAPACK_FUNC_UNDERSCORES))
    print_error("Both SUNDIALS_LAPACK_FUNC_CASE and SUNDIALS_LAPACK_FUNC_UNDERSCORES must be set.")
  endif()

  string(TOUPPER ${SUNDIALS_LAPACK_FUNC_CASE} _case)
  string(TOUPPER ${SUNDIALS_LAPACK_FUNC_UNDERSCORES} _underscores)

  # Set the C preprocessor macro
  set(LAPACK_MANGLE_MACRO "#define SUNDIALS_LAPACK_FUNC(name,NAME)")

  if(_case MATCHES "LOWER")
    set(LAPACK_MANGLE_MACRO "${LAPACK_MANGLE_MACRO} name")
  elseif(_case MATCHES "UPPER")
    set(LAPACK_MANGLE_MACRO "${LAPACK_MANGLE_MACRO} NAME")
  else()
    print_error("Invalid SUNDIALS_LAPACK_FUNC_CASE option.")
  endif()

  if(SUNDIALS_LAPACK_FUNC_SUFFIX)
    set(LAPACK_MANGLE_MACRO "${LAPACK_MANGLE_MACRO} ## ${SUNDIALS_LAPACK_FUNC_SUFFIX}")
  endif()

  if(_underscores MATCHES "ONE")
    set(LAPACK_MANGLE_MACRO "${LAPACK_MANGLE_MACRO} ## _")
  elseif(_underscores MATCHES "TWO")
    set(LAPACK_MANGLE_MACRO "${LAPACK_MANGLE_MACRO} ## __")
  elseif(NOT (_underscores MATCHES "NONE"))
    print_error("Invalid SUNDIALS_LAPACK_FUNC_UNDERSCORES option.")
  endif()

  # Test timers with a simple program
  set(LAPACK_TEST_DIR ${PROJECT_BINARY_DIR}/LapackTest)
  file(MAKE_DIRECTORY ${LAPACK_TEST_DIR})

  # Create a CMakeLists.txt file which will generate the test executable
  file(WRITE ${LAPACK_TEST_DIR}/CMakeLists.txt
    "cmake_minimum_required(VERSION ${CMAKE_VERSION})\n"
    "project(ltest C)\n"
    "set(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "set(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "set(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
    "set(CMAKE_C_STANDARD ${CMAKE_C_STANDARD})\n"
    "set(CMAKE_C_EXTENSIONS ${CMAKE_C_EXTENSIONS})\n"
    "set(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "set(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "set(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "set(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "set(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "add_executable(ltest ltest.c)\n"
    "target_link_libraries(ltest \"${LAPACK_LIBRARIES}\")\n")

  # Create a simple C source for testing
  file(WRITE ${LAPACK_TEST_DIR}/ltest.c
    "${LAPACK_MANGLE_MACRO}\n"
    "#define dcopy_mangled SUNDIALS_LAPACK_FUNC(dcopy, DCOPY)\n"
    "#define dgetrf_mangled SUNDIALS_LAPACK_FUNC(dgetrf, DGETRF)\n"
    "extern void dcopy_mangled(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);\n"
    "extern void dgetrf_magnled(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);\n"
    "int main(){\n"
    "int n=1;\n"
    "double x, y;\n"
    "dcopy_mangled(&n, &x, &n, &y, &n);\n"
    "dgetrf_mangled(&n, &n, &x, &n, &n, &n);\n"
    "return 0;\n"
    "}\n")

  # To ensure we do not use stuff from the previous attempts,
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${LAPACK_TEST_DIR}/CMakeFiles)

  # Use TRY_COMPILE to make the target
  try_compile(COMPILE_OK ${LAPACK_TEST_DIR} ${LAPACK_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

endmacro()


if(NOT LAPACK_WORKS)

  # Test the current settings first
  test_lapack_name_mangling()

  # Test the possible options
  if(NOT COMPILE_OK)

    foreach(case "LOWER" "UPPER")
      foreach(underscores "NONE" "ONE" "TWO")

        # Overwrite the cache variable values
        set(SUNDIALS_LAPACK_FUNC_CASE "${case}" CACHE STRING
          "Case of LAPACK function names" FORCE)

        set(SUNDIALS_LAPACK_FUNC_UNDERSCORES "${underscores}" CACHE STRING
          "Number of underscores appended to LAPACK function names (none/one/two)"
          FORCE)

        test_lapack_name_mangling()
        if(COMPILE_OK)
          break()
        endif()

      endforeach()
    endforeach()

  endif()

  # Process test result
  if(COMPILE_OK)

    message(STATUS "Checking if LAPACK works with SUNDIALS... OK")
    set(LAPACK_WORKS TRUE CACHE BOOL "LAPACK works with SUNDIALS as configured" FORCE)

    # get path to LAPACK library to use in generated makefiles for examples, if
    # LAPACK_LIBRARIES contains multiple items only use the path of the first entry
    list(LENGTH LAPACK_LIBRARIES len)
    if(len EQUAL 1)
      get_filename_component(LAPACK_LIBRARY_DIR ${LAPACK_LIBRARIES} PATH)
    else()
      list(GET LAPACK_LIBRARIES 0 TMP_LAPACK_LIBRARIES)
      get_filename_component(LAPACK_LIBRARY_DIR ${TMP_LAPACK_LIBRARIES} PATH)
    endif()

  else()

    message(STATUS "Checking if LAPACK works with SUNDIALS... FAILED")
    set(LAPACK_WORKS FALSE CACHE BOOL "LAPACK does not work with SUNDIALS as configured" FORCE)
    print_error("SUNDIALS interface to LAPACK is not functional.")

  endif()

else()

  message(STATUS "Skipped LAPACK tests, assuming LAPACK works with SUNDIALS. Set LAPACK_WORKS=FALSE to (re)run compatibility test.")

endif()
