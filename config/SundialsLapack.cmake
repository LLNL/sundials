# ---------------------------------------------------------------
# $Revision: 1.1 $
# $Date: 2009-02-17 02:58:46 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2008, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# BLAS/LAPACK tests for SUNDIALS CMake-based configuration.
#
# 

SET(LAPACK_FOUND FALSE)

# If LAPACK libraries are undefined, try to find them (if we have
# a working Fortran compiler) or look for them in the most
# obvious place...
if(NOT LAPACK_LIBRARIES)
  if(F77_FOUND)
    include(FindLAPACK)
  else(F77_FOUND)
    find_library(LAPACK_LIBRARIES
      NAMES lapack
      PATHS /usr/lib /usr/local/lib
      "$ENV{ProgramFiles}/LAPACK/Lib"
      )
  endif(F77_FOUND)
endif(NOT LAPACK_LIBRARIES)
# If using a GNU C compiler, it is quite likely we'll want LAPACK_LINKER_FLAGS
# to include -lg2c (if not already present)
if(CMAKE_COMPILER_IS_GNUCC AND NOT LAPACK_LINKER_FLAGS MATCHES "g2c")
  set(LAPACK_LINKER_FLAGS "${LAPACK_LINKER_FLAGS} -lg2c")
endif(CMAKE_COMPILER_IS_GNUCC AND NOT LAPACK_LINKER_FLAGS MATCHES "g2c")
# If we have the LAPACK libraries, test them
if(LAPACK_LIBRARIES)
  message(STATUS "Looking for LAPACK libraries... OK")
  # Create the LapackTest directory
  set(LapackTest_DIR ${PROJECT_BINARY_DIR}/LapackTest)
  file(MAKE_DIRECTORY ${LapackTest_DIR})
  # Create a CMakeLists.txt file 
  file(WRITE ${LapackTest_DIR}/CMakeLists.txt
    "CMAKE_MINIMUM_REQUIRED(VERSION 2.2)\n"
    "PROJECT(ltest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
    "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n"
    "SET(CMAKE_C_FLAGS_RELEASE \"${CMAKE_C_FLAGS_RELEASE}\")\n"
    "SET(CMAKE_C_FLAGS_DEBUG \"${CMAKE_C_FLAGS_DEBUG}\")\n"
    "SET(CMAKE_C_FLAGS_RELWITHDEBUGINFO \"${CMAKE_C_FLAGS_RELWITHDEBUGINFO}\")\n"
    "SET(CMAKE_C_FLAGS_MINSIZE \"${CMAKE_C_FLAGS_MINSIZE}\")\n"
    "SET(CMAKE_EXE_LINKER_FLAGS \"\${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${LAPACK_LIBRARIES})\n")    
  # Create a C source file which calls a Blas function (dcopy) and an Lapack function (dgetrf)
  file(WRITE ${LapackTest_DIR}/ltest.c
    "${F77_MANGLE_MACRO1}\n"
    "#define dcopy_f77 SUNDIALS_F77_FUNC(dcopy, DCOPY)\n"
    "#define dgetrf_f77 SUNDIALS_F77_FUNC(dgetrf, DGETRF)\n"
    "extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);\n"
    "extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);\n"
    "int main(){\n"
    "int n=1;\n"
    "double x, y;\n"
    "dcopy_f77(&n, &x, &n, &y, &n);\n"
    "dgetrf_f77(&n, &n, &x, &n, &n, &n);\n"
    "return(0);\n"
    "}\n")
  # Attempt to link the "ltest" executable
  try_compile(LTEST_OK ${LapackTest_DIR} ${LapackTest_DIR}
    ltest OUTPUT_VARIABLE MY_OUTPUT)    
  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  file(REMOVE_RECURSE ${LapackTest_DIR}/CMakeFiles)
  # Process test result
  if(LTEST_OK)
    message(STATUS "Checking if Lapack works... OK")
    set(LAPACK_FOUND TRUE)
  else(LTEST_OK)
    message(STATUS "Checking if Lapack works... FAILED")
  endif(LTEST_OK)
else(LAPACK_LIBRARIES)
  message(STATUS "Looking for LAPACK libraries... FAILED")
endif(LAPACK_LIBRARIES)
