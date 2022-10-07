# ---------------------------------------------------------------
# Programmer(s): Radu Serban, David Gardner, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ---------------------------------------------------------------
# Module which enables Fortran and tests for support of necessary
# compiler features for the current SUNDIALS configuration.
# Will define the variables:
#   Fortran_FOUND - TRUE if a Fortran compiler is found
#   F2003_FOUND   - TRUE if the Fortran compiler supports the
#                   Fortran 2003 standard
# ---------------------------------------------------------------

# If the Fortran compiler flags are set using environemnt variables (i.e.,
# CMAKE_Fortran_FLAGS is not set), then check if both FFLAGS and FCFLAGS are
# set. If both are set and not the same then a fatal error occurs.
#
# NOTE: This check must occur before 'enable_language(Fortran)' as it will use
# the value of FFLAGS to set CMAKE_Fortran_FLAGS
set(ENV_FFLAGS "$ENV{FFLAGS}")
set(ENV_FCFLAGS "$ENV{FCFLAGS}")

# check if environment variables are used and CMAKE_Fortran_FLAGS is not
if ((NOT "${ENV_FFLAGS}" STREQUAL "") AND (NOT "${ENV_FCFLAGS}" STREQUAL "")
    AND ("${CMAKE_Fortran_FLAGS}" STREQUAL ""))

  # check if environment variables are equal
  if (NOT "${ENV_FFLAGS}" STREQUAL "${ENV_FCFLAGS}")
    print_error("FFLAGS='${ENV_FFLAGS}' and FCFLAGS='${ENV_FCFLAGS}' are both set but are not equal.")
  endif()
endif()

# -----------------------------------------------------------------------------
# Enable Fortran
# -----------------------------------------------------------------------------
enable_language(Fortran)
set(Fortran_FOUND TRUE)

# Enable preprocessing Fortran code. With older versions of CMake is this
# handled in SundialsAddLibrary.cmake by adding a compiler option.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
  set(CMAKE_Fortran_PREPROCESS ON)
endif()

# -----------------------------------------------------------------------------
# Check if Fortran 2003 is supported
# -----------------------------------------------------------------------------
if(BUILD_FORTRAN_MODULE_INTERFACE)
  if(NOT F2003_FOUND)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003")

    set(F2003Test_DIR ${PROJECT_BINARY_DIR}/F2003Test_DIR)
    file(MAKE_DIRECTORY ${F2003Test_DIR})

    # Create a CMakeLists.txt file
    file(WRITE ${F2003Test_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ftest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
      "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ftest ftest.f90)\n")

    # Create a Fortran source file which tries to use iso_c_binding
    file(WRITE ${F2003Test_DIR}/ftest.f90
      "program main\n"
      "use, intrinsic :: iso_c_binding\n"
      "end program main\n")

    # Attempt compile the executable
    try_compile(FTEST_OK ${F2003Test_DIR} ${F2003Test_DIR}
      ftest OUTPUT_VARIABLE COMPILE_OUTPUT)

    # To ensure we do not use stuff from the previous attempts,
    # we must remove the CMakeFiles directory.
    file(REMOVE_RECURSE ${F2003Test_DIR}/CMakeFiles)

    if(FTEST_OK)
      message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003 -- yes")
      set(F2003_FOUND TRUE CACHE BOOL "${CMAKE_Fortran_COMPILER} supports F2003" FORCE)
    else()
      message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003 -- no")
      message(STATUS "Check output:")
      message("${COMPILE_OUTPUT}")
      print_error("BUILD_FORTRAN_MODULE_INTERFACE is set to ON, but the CMAKE_Fortran_COMPILER does not support F2003")
    endif()
  else()
    message(STATUS "Skipped F2003 tests, assuming ${CMAKE_Fortran_COMPILER} supports the f2003 standard. To rerun the F2003 tests, set F2003_FOUND to FALSE.")
  endif()
endif()

# Put all F2003 modules into one build directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran")

# ---------------------------------------------------------------
# Determining the name-mangling scheme if needed
# ---------------------------------------------------------------
# In general, names of symbols with and without underscore may be mangled
# differently (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__),
# we have to consider both cases.
#
# Method:
#  1) create a library from a Fortran source file which defines a function "mysub"
#  2) attempt to link with this library a C source file which calls the "mysub"
#     function using various possible schemes (6 different schemes, corresponding
#     to all combinations lower/upper case and none/one/two underscores).
#  3) define the name-mangling scheme based on the test that was successful.
#
# On exit, if we were able to infer the scheme, the variables
# CMAKE_Fortran_SCHEME_NO_UNDERSCORES and CMAKE_Fortran_SCHEME_WITH_UNDERSCORES
# contain the mangled names for "mysub" and "my_sub", respectively.
# ---------------------------------------------------------------
if(NEED_FORTRAN_NAME_MANGLING)

  set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "")
  set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "")

  # Create the FortranTest directory
  set(FortranTest_DIR ${PROJECT_BINARY_DIR}/FortranTest)
  file(MAKE_DIRECTORY ${FortranTest_DIR})

  # Create a CMakeLists.txt file which will generate the "flib" library
  # and an executable "ftest"
  file(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
    "PROJECT(ftest Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
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

    # Infer Fortran name-mangling scheme for symbols WITHOUT underscores.
    # Overwrite CMakeLists.txt with one which will generate the "ctest1" executable
    file(WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ctest1 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest1 ctest1.c)\n"
      "FIND_LIBRARY(FLIB flib \"${FortranTest_DIR}\")\n"
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
      file(WRITE ${FortranTest_DIR}/ctest1.c
        "extern void ${opt}();\n"
        "int main(){${opt}();return(0);}\n")
      # Use TRY_COMPILE to make the "ctest1" executable from the current C source
      # and linking to the previously created "flib" library.
      try_compile(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
        ctest1 OUTPUT_VARIABLE MY_OUTPUT)
      # Write output compiling the test code
      file(WRITE ${FortranTest_DIR}/ctest1_${opt}.out "${MY_OUTPUT}")
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
      endif()
    endwhile(${iopt} LESS ${imax})

    # Infer Fortran name-mangling scheme for symbols WITH underscores.
    # Practically a duplicate of the previous steps.
    file(WRITE ${FortranTest_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ctest2 C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_C_COMPILER \"${CMAKE_C_COMPILER}\")\n"
      "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
      "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ctest2 ctest2.c)\n"
      "FIND_LIBRARY(FLIB flib \"${FortranTest_DIR}\")\n"
      "TARGET_LINK_LIBRARIES(ctest2 \${FLIB})\n")

    set(options my_sub my_sub_ my_sub__ MY_SUB MY_SUB_ MY_SUB__)
    list(LENGTH options imax)
    set(iopt 0)
    while(${iopt} LESS ${imax})
      list(GET options ${iopt} opt)
      file(WRITE ${FortranTest_DIR}/ctest2.c
        "extern void ${opt}();\n"
        "int main(){${opt}();return(0);}\n")
      try_compile(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
        ctest2 OUTPUT_VARIABLE MY_OUTPUT)
      file(WRITE ${FortranTest_DIR}/ctest2_${opt}.out "${MY_OUTPUT}")
      file(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
      if(CTEST_OK)
        set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES ${opt})
        set(iopt ${imax})
      else(CTEST_OK)
        math(EXPR iopt ${iopt}+1)
      endif()
    endwhile(${iopt} LESS ${imax})

    # If a name-mangling scheme was found set the C preprocessor macros to use
    # that scheme. Otherwise default to lower case with one underscore.
    if(CMAKE_Fortran_SCHEME_NO_UNDERSCORES AND CMAKE_Fortran_SCHEME_WITH_UNDERSCORES)
      message(STATUS "Determining Fortran name-mangling scheme... OK")
    else()
      message(STATUS "Determining Fortran name-mangling scheme... DEFAULT")
      set(CMAKE_Fortran_SCHEME_NO_UNDERSCORES "mysub_")
      set(CMAKE_Fortran_SCHEME_WITH_UNDERSCORES "my_sub_")
    endif()

    # Symbols NO underscores
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## __")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
      set(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## __")
    endif()

    # Symbols WITH underscores
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub_")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "my_sub__")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## __")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB_")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## _")
    endif()
    if(${CMAKE_Fortran_SCHEME_WITH_UNDERSCORES} MATCHES "MY_SUB__")
      set(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## __")
    endif()

    # name-mangling scheme has been set
    set(NEED_FORTRAN_NAME_MANGLING FALSE)
  else(FTEST_OK)
    message(STATUS "Determining Fortran name-mangling scheme... FAILED")
  endif()

endif()
