# ---------------------------------------------------------------
# $Revision:  $
# $Date:  $
# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# PETSc tests for SUNDIALS CMake-based configuration.
#    - loosely based on SundialsLapack.cmake
# 

SET(PETSC_FOUND FALSE)

# set PETSC_LIBRARIES
include(FindPETSc)
# If we have the PETSC libraries, test them
if(PETSC_LIBRARY)
  message(STATUS "Looking for PETSc libraries...")
  # Create the PETSCTest directory
  set(PETSCTest_DIR ${PROJECT_BINARY_DIR}/PETSCTest)
  file(MAKE_DIRECTORY ${PETSCTest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${PETSCTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${PETSC_LIBRARY})\n")    
# Create a C source file which calls a PETSC function
  file(WRITE ${PETSCTest_DIR}/ltest.c
    "\#include \"petscvec.h\"\n"
    "int main(){\n"
    "Vec x;\n"
    "VecCreate(PETSC_COMM_WORLD, &x);\n" 
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${PETSCTest_DIR} ${PETSCTest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${PETSCTest_DIR}/CMakeFiles)
  # Process test result
#PRINT_WARNING("LTEST_OK" "${LTEST_OK}")
  if(LTEST_OK)
#PRINT_WARNING("x SundialsPETSC.cmake PETSC_LIBRARIES" "${PETSC_LIBRARIES}")
    message(STATUS "Checking if PETSc works... OK")
    set(PETSC_FOUND TRUE)
    #print_warning("PETSC_FOUND" "${PETSC_FOUND}")
  else(LTEST_OK)
    message(STATUS "Checking if PETSc works... FAILED")
  endif(LTEST_OK)
else(PETSC_LIBRARY)
#PRINT_WARNING("y SundialsPETSC.cmake PETSC_LIBRARIES" "${PETSC_LIBRARIES}")
  message(STATUS "Looking for PETSc libraries... FAILED")
endif(PETSC_LIBRARY)
