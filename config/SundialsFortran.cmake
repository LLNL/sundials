# ---------------------------------------------------------------
# $Revision: 4511 $
# $Date: 2015-06-29 16:48:09 -0700 (Mon, 29 Jun 2015) $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2008, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Fortran-related tests for SUNDIALS CMake-based configuration.
#
# Determining the name-mangling scheme
# ------------------------------------
# In general, names of symbols with and without underscore may be mangled 
# differently (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__),
# we have to consider both cases.
# Method:
#  1) create a library from a Fortran source file which defines a function "mysub"
#  2) attempt to link with this library a C source file which calls the "mysub"
#     function using various possible schemes (6 different schemes, corresponding
#     to all combinations lower/upper case and none/one/two underscores)
#  3) define the name-mangling scheme based on the test that was successful.
# On exit, if we were able to infer the scheme, the variables
# CMAKE_Fortran_SCHEME_NO_UNDERSCORES and CMAKE_Fortran_SCHEME_WITH_UNDERSCORES
# contain the mangled names for "mysub" and "my_sub", respectively.

set(F77_FOUND FALSE)
set(F77SCHEME_FOUND FALSE)

set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "")
set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "")

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
  # Create a CMakeLists.txt file which will generate the "flib" library
  # and an executable "ftest"
  file(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
    "PROJECT(ftest Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
    "ADD_LIBRARY(flib flib.f)\n"
    "ADD_EXECUTABLE(ftest ftest.f)\n"
    "TARGET_LINK_LIBRARIES(ftest flib)\n")
  # Create the Fortran source flib.f which defines two subroutines, "mysub" and "my_sub"
  file(WRITE ${FortranTest_DIR}/flib.f
    "        SUBROUTINE mysub\n"
    "        RETURN\n"
    "        END\n"
    "        SUBROUTINE my_sub\n"
    "        RETURN\n"
    "        END\n")
  # Create the Fortran source ftest.f which calls "mysub" and "my_sub"
  file(WRITE ${FortranTest_DIR}/ftest.f
    "        PROGRAM ftest\n"
    "        CALL mysub()\n"
    "        CALL my_sub()\n"
    "        END\n")
  # Use TRY_COMPILE to make the targets "flib" and "ftest"
  try_compile(FTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
    ftest OUTPUT_VARIABLE MY_OUTPUT)
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
  # Proceed based on test results
  if(FTEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran program... OK")
    set(F77_FOUND TRUE)
    # Infer Fortran name-mangling scheme for symbols WITHOUT underscores.
    # Overwrite CMakeLists.txt with one which will generate the "ctest1" executable
    file(WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(ctest1 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest1 ctest1.c)\n"
      "FIND_LIBRARY(FLIB flib ${FortranTest_DIR})\n"
      "TARGET_LINK_LIBRARIES(ctest1 \${FLIB})\n")
    # Define the list "options" of all possible schemes that we want to consider
    # Get its length and initialize the counter "iopt" to zero
    set(options mysub mysub_ mysub__ MYSUB MYSUB_ MYSUB__)
    list(LENGTH options imax)
    set(iopt 0)
    # We will attempt to sucessfully generate the "ctest1" executable as long as
    # there still are entries in the "options" list
    while(${iopt} LESS ${imax})   
      # Get the current list entry (current scheme)
      list(GET options ${iopt} opt)
      # Generate C source which calls the "mysub" function using the current scheme
      file(WRITE ${FortranTest_DIR}/ctest1.c "int main(){${opt}();return(0);}\n")
      # Use TRY_COMPILE to make the "ctest1" executable from the current C source
      # and linking to the previously created "flib" library.
      try_compile(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
        ctest1 OUTPUT_VARIABLE MY_OUTPUT)
      # To ensure we do not use stuff from the previous attempts, 
      # we must remove the CMakeFiles directory.
      file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
      # Test if we successfully created the "ctest" executable.
      # If yes, save the current scheme, and set the counter "iopt" to "imax" 
      # so that we exit the while loop.
      # Otherwise, increment the counter "iopt" and go back in the while loop.
      if(CTEST_OK)
        set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES ${opt})
        set(iopt ${imax})
      else(CTEST_OK)
        math(EXPR iopt ${iopt}+1)
      endif(CTEST_OK)
    endwhile(${iopt} LESS ${imax})   
    # Infer Fortran name-mangling scheme for symbols WITH underscores.
    # Practically a duplicate of the previous steps.
    file(WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION 2.4)\n"
      "PROJECT(ctest2 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest2 ctest2.c)\n"
      "FIND_LIBRARY(FLIB flib ${FortranTest_DIR})\n"
      "TARGET_LINK_LIBRARIES(ctest2 \${FLIB})\n")
    set(options my_sub my_sub_ my_sub__ MY_SUB MY_SUB_ MY_SUB__)
    list(LENGTH options imax)
    set(iopt 0)
    while(${iopt} LESS ${imax})   
      list(GET options ${iopt} opt)
      file(WRITE ${FortranTest_DIR}/ctest2.c "int main(){${opt}();return(0);}\n")
      try_compile(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
        ctest2 OUTPUT_VARIABLE MY_OUTPUT)
      file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
      if(CTEST_OK)
        set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES ${opt})
        set(iopt ${imax})
      else(CTEST_OK)
        math(EXPR iopt ${iopt}+1)
      endif(CTEST_OK)
    endwhile(${iopt} LESS ${imax})   
    # Proceed based on whether the previous tests were successfull or not
    if(CMAKE_Fortran_SCHEME_NO_UNDERSCORES AND CMAKE_Fortran_SCHEME_WITH_UNDERSCORES)
      message(STATUS "Determining Fortran name-mangling scheme... OK")
      set(F77SCHEME_FOUND TRUE)
    else(CMAKE_Fortran_SCHEME_NO_UNDERSCORES AND CMAKE_Fortran_SCHEME_WITH_UNDERSCORES)
      message(STATUS "Determining Fortran name-mangling scheme... FAILED")
    endif(CMAKE_Fortran_SCHEME_NO_UNDERSCORES AND CMAKE_Fortran_SCHEME_WITH_UNDERSCORES)
  else(FTEST_OK)
    message(STATUS "Trying to compile and link a simple Fortran program... FAILED")
  endif(FTEST_OK)
else(CMAKE_Fortran_COMPILER)
  message(STATUS "Searching for a Fortran compiler... FAILED")
endif(CMAKE_Fortran_COMPILER)

