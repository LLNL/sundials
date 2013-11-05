# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2013-02-14 $
# ---------------------------------------------------------------
# Programmer:  Steven Smith @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# KLU tests for SUNDIALS CMake-based configuration.
#    - loosely based on SundialsLapack.cmake
# 

SET(KLU_FOUND FALSE)

# set KLU_LIBRARIES
include(FindKLU)
# If we have the KLU libraries, test them
if(KLU_LIBRARIES)
  message(STATUS "Looking for KLU libraries...")
  # Create the KLUTest directory
  set(KLUTest_DIR ${PROJECT_BINARY_DIR}/KLUTest)
  file(MAKE_DIRECTORY ${KLUTest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${KLUTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.2)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "INCLUDE_DIRECTORIES(${KLU_INCLUDE_DIR})\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${KLU_LIBRARIES})\n")    
# Create a C source file which calls a KLU function
# SGS TODO what is a simple KLU method to invoke?
  file(WRITE ${KLUTest_DIR}/ltest.c
    "\#include \"klu.h\"\n"
    "int main(){\n"
    "klu_common Common;\n"
    "klu_defaults (&Common);\n" 
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${KLUTest_DIR} ${KLUTest_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${KLUTest_DIR}/CMakeFiles)
  # Process test result
#PRINT_WARNING("LTEST_OK" "${LTEST_OK}")
  if(LTEST_OK)
#PRINT_WARNING("x SundialsKLU.cmake KLU_LIBRARIES" "${KLU_LIBRARIES}")
    message(STATUS "Checking if KLU works... OK")
    set(KLU_FOUND TRUE)
    #print_warning("KLU_FOUND" "${KLU_FOUND}")
  else(LTEST_OK)
    message(STATUS "Checking if KLU works... FAILED")
  endif(LTEST_OK)
else(KLU_LIBRARIES)
#PRINT_WARNING("y SundialsKLU.cmake KLU_LIBRARIES" "${KLU_LIBRARIES}")
  message(STATUS "Looking for KLU libraries... FAILED")
endif(KLU_LIBRARIES)
