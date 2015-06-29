# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
# ---------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Fortran90-related tests for SUNDIALS CMake-based configuration.

set(F90_FOUND FALSE)
include(CMakeDetermineFortranCompiler)

if(CMAKE_Fortran_COMPILER)
  message(STATUS "Searching for a Fortran compiler... ${CMAKE_Fortran_COMPILER}")
  # Enable the language for next steps
  enable_language(Fortran)
  mark_as_advanced(CLEAR 
    CMAKE_Fortran_COMPILER
    CMAKE_Fortran_FLAGS
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_MINSIZEREL
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_RELWITHDEB)
  # Create the FortranTest directory
  set(FortranTest_DIR ${PROJECT_BINARY_DIR}/FortranTest)
  file(MAKE_DIRECTORY ${FortranTest_DIR})
  # Create a CMakeLists.txt file which will generate the "f90lib" library
  # and an executable "f90test"
  file(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(f90test Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
    "ADD_LIBRARY(f90lib f90lib.f90)\n"
    "ADD_EXECUTABLE(f90test f90test.f90)\n"
    "TARGET_LINK_LIBRARIES(f90test f90lib)\n")
  # Create the Fortran source f90lib.f90 which defines two subroutines, "mysub" and "my_sub"
  file(WRITE ${FortranTest_DIR}/f90lib.f90
    "subroutine mysub\n"
    "  return\n"
    "end\n"
    "subroutine my_sub\n"
    "  return\n"
    "end\n")
  # Create the Fortran source f90test.f90 which calls "mysub" and "my_sub"
  file(WRITE ${FortranTest_DIR}/f90test.f90
    "program f90test\n"
    "  call mysub()\n"
    "  call my_sub()\n"
    "end\n")
  # Use TRY_COMPILE to make the targets "f90lib" and "f90test"
  try_compile(F90TEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
    f90test OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
  # Proceed based on test results
  if(F90TEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran90 program... OK")
    set(F90_FOUND TRUE)
  else(F90TEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran90 program... FAILED")
  endif(F90TEST_OK)
else(CMAKE_Fortran_COMPILER)
  message(STATUS "Searching for a Fortran compiler... FAILED")
endif(CMAKE_Fortran_COMPILER)

