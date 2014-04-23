# ---------------------------------------------------------------
# $Revision$
# $Date$
# ---------------------------------------------------------------
# Programmer:  Eddy Banks @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2013, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# SUPERLUMT tests for SUNDIALS CMake-based configuration.
#    - loosely based on SundialsLapack.cmake
# 

#print_warning("SundialsSUPERLUMT.cmake 1 SUPERLUMT_FOUND" "${SUPERLUMT_FOUND}")
SET(SUPERLUMT_FOUND FALSE)
#print_warning("SundialsSUPERLUMT.cmake 2 SUPERLUMT_FOUND" "${SUPERLUMT_FOUND}")

# set SUPERLUMT_LIBRARIES
include(FindSUPERLUMT)
# If we have the SUPERLUMT libraries, test them
#print_warning("SundialsSUPERLUMT.cmake 3: SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
#print_warning("SundialsSUPERLUMT.cmake 4: SUPERLUMT_BLAS_LIBRARIES" "${SUPERLUMT_BLAS_LIBRARIES}")
if(SUPERLUMT_LIBRARIES)
  #print_warning("SundialsSUPERLUMT.cmake 5 SUPERLUMT_FOUND" "${SUPERLUMT_FOUND}")
  #PRINT_WARNING("SundialsSUPERLUMT.cmake 6 SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
  message(STATUS "Looking for SUPERLUMT libraries... OK")
  # Create the SUPERLUMT_TEST directory
  set(SUPERLUMT_TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUMT_TEST)
  file(MAKE_DIRECTORY ${SUPERLUMT_TEST_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${SUPERLUMT_TEST_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.2)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${SUPERLUMT_LIBRARIES})\n")    
# TODO: Eddy - fix this test
# Create a C source file which calls a SUPERLUMT function
  file(WRITE ${SUPERLUMT_TEST_DIR}/ltest.c
    "int main(){\n"
    "int n=1;\n"
    "double x, y;\n"
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${SUPERLUMT_TEST_DIR} ${SUPERLUMT_TEST_DIR} ltest OUTPUT_VARIABLE MY_OUTPUT)
      
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${SUPERLUMT_TEST_DIR}/CMakeFiles)
  # Process test result
#PRINT_WARNING("LTEST_OK" "${LTEST_OK}")
  if(LTEST_OK)
#PRINT_WARNING("x SundialsSUPERLUMT.cmake SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
    message(STATUS "Checking if SUPERLUMT works... OK")
    set(SUPERLUMT_FOUND TRUE)
    #print_warning("SUPERLUMT_FOUND" "${SUPERLUMT_FOUND}")
  else(LTEST_OK)
    message(STATUS "Checking if SUPERLUMT works... FAILED")
  endif(LTEST_OK)
else(SUPERLUMT_LIBRARIES)
#PRINT_WARNING("y SundialsSUPERLUMT.cmake SUPERLUMT_LIBRARIES" "${SUPERLUMT_LIBRARIES}")
  message(STATUS "Looking for SUPERLUMT libraries... FAILED")
endif(SUPERLUMT_LIBRARIES)
