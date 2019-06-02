# ------------------------------------------------------------------------------
# Programmer:  Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Sundials module to find and configure SuperLU_DIST correctly.
# -----------------------------------------------------------------------------

### This is only set if running GUI - simply return first time enabled
if(SUPERLUDIST_DISABLED)
  set(SUPERLUDIST_DISABLED FALSE CACHE INTERNAL "GUI - SUPERLUDIST now enabled" FORCE)
  return()
endif()

# --- Check for dependencies --- #

# SuperLU_DIST only supports double precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  PRINT_ERROR("SuperLU_DIST is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# SuperLU_DIST requires MPI, so make sure it was found.
if(NOT (MPI_C_FOUND AND MPI_CXX_FOUND))
  PRINT_ERROR("SuperLU_DIST requires MPI but it was not found.")
endif()

# SuperLU_DIST OpenMP node parallelism is on, make sure OpenMP as found and is
# at least version 4.5
if(SUPERLUDIST_OpenMP)

  if(NOT OPENMP_FOUND)
    PRINT_ERROR("SUPERLUDIST_OpenMP is set to ON but OpenMP was not found.")
  elseif(NOT OPENMP45_FOUND)
    string(CONCAT ERRSTR "SuperLUDIST requires OpenMP 4.5+ but it was not found. "
      "Either use CMake 3.9+, or if you are sure OpenMP 4.5+ is available "
      "set the SKIP_OPENMP_DEVICE_CHECK advanced option to ON.")
    PRINT_ERROR(${ERRSTR})
  endif()

endif()

# --- Find SuperLU_DIST and test it --- #

# Try and find SuperLU_DIST now
find_package(SUPERLUDIST REQUIRED)

# If we have the SuperLUDIST libraries, test them
if(SUPERLUDIST_FOUND)
  message(STATUS "Checking if SuperLUDIST works... ")

  # Check index size
  if(NOT (SUNDIALS_INDEX_SIZE STREQUAL SUPERLUDIST_INDEX_SIZE))
    set(_err_msg_string "SuperLU_DIST not functional due to index size mismatch:\n")
    string(APPEND _err_msg_string "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, but SuperLU_DIST was built with ${SUPERLUDIST_INDEX_SIZE}-bit indices\n")
    string(APPEND _err_msg_string "SUPERLUDIST_INCLUDE_DIR: ${SUPERLUDIST_INCLUDE_DIR}\n")
    PRINT_ERROR("${_err_msg_string}")
  endif()

  # Create the SUPERLUDIST_TEST directory
  set(SUPERLUDIST_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUDIST_Test)
  file(MAKE_DIRECTORY ${SUPERLUDIST_TEST_DIR})

  # Create a CMakeLists.txt file
  file(WRITE ${SUPERLUDIST_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 3.1.3)\n"
    "PROJECT(ltest CXX)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_CXX_COMPILER \"${MPI_CXX_COMPILER}\")\n"
    "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${SUPERLUDIST_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.cpp)\n"
    "LINK_LIBRARIES(${SUPERLUDIST_LIBRARIES})\n")

  # Create a C source file which calls a SuperLUDIST function
  # and also prints the size of the indices used.
  file(WRITE ${SUPERLUDIST_TEST_DIR}/ltest.cpp
    "\#include <superlu_ddefs.h>\n"
    "int main(){\n"
    "SuperMatrix *A;\n"
    "NRformat_loc *Astore;\n"
    "A = NULL;\n"
    "Astore = NULL;\n"
    "if (A != NULL || Astore != NULL) return(1);\n"
    "else return(0);\n"
    "}\n")

  try_compile(COMPILE_OK ${SUPERLUDIST_TEST_DIR} ${SUPERLUDIST_TEST_DIR} ltest
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if SuperLUDIST works... OK")
  else()
    message(STATUS "Checking if SuperLUDIST works... FAILED")
    if(NOT COMPILE_OK)
      message(STATUS "Check output: ")
      message("${COMPILE_OUTPUT}")
      PRINT_ERROR("SuperLUDIST not functional - support will not be provided.")
    endif()
  endif()

  # sundials_config.h symbols
  set(SUNDIALS_SUPERLUDIST TRUE)

else()

  # sundials_config.h symbols
  set(SUNDIALS_SUPERLUDIST FALSE)

endif()
