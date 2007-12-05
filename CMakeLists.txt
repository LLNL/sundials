# ---------------------------------------------------------------
# $Revision: 1.5 $
# $Date: 2007-12-05 23:33:24 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Top level CMakeLists.txt for SUNDIALS (for cmake build system)


# -------------------------------------------------------------
# Initial commands
# -------------------------------------------------------------

# Require a fairly recent cmake version

CMAKE_MINIMUM_REQUIRED(VERSION 2.2)

# Project SUNDIALS (initially only C supported)

PROJECT(sundials C)

# Set some variables with info on the SUNDIALS project

SET(PACKAGE_BUGREPORT "radu@llnl.gov")
SET(PACKAGE_NAME "SUNDIALS")
SET(PACKAGE_STRING "SUNDIALS 2.4.0")
SET(PACKAGE_TARNAME "sundials")
SET(PACKAGE_VERSION "2.4.0")

# Prohibit in-source build

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Hide/show some cache variables

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CLEAR
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEB)

# Specify the VERSION and SOVERSION for shared libraries

SET(cvodelib_VERSION "1.0.0")
SET(cvodelib_SOVERSION "1")

SET(cvodeslib_VERSION "2.0.0")
SET(cvodeslib_SOVERSION "2")

SET(idalib_VERSION "2.0.0")
SET(idalib_SOVERSION "2")

SET(idaslib_VERSION "0.0.0")
SET(idaslib_SOVERSION "0")

SET(kinsollib_VERSION "1.0.0")
SET(kinsollib_SOVERSION "1")

SET(cpodeslib_VERSION "0.0.0")
SET(cpodeslib_SOVERSION "0")

SET(nveclib_VERSION "0.0.2")
SET(nveclib_SOVERSION "0")

# Specify the location of additional CMAKE modules

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# -------------------------------------------------------------
# MACRO definitions
# -------------------------------------------------------------

# Macros to hide/show cached variables.
# These two macros can be used to "hide" or "show" in the
# list of cached variables various variables and/or options 
# that depend on other options.
# Note that once a variable is modified, it will preserve its
# value (hidding it merely makes it internal)

MACRO(HIDE_VARIABLE var)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE INTERNAL "")
  ENDIF(DEFINED ${var})
ENDMACRO(HIDE_VARIABLE)

MACRO(SHOW_VARIABLE var type doc default)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE "${type}" "${doc}" FORCE)
  ELSE(DEFINED ${var})
    SET(${var} "${default}" CACHE "${type}" "${doc}")
  ENDIF(DEFINED ${var})
ENDMACRO(SHOW_VARIABLE)

# Macros to append a common suffix or prefix to the elements of a list

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

# Macro to print warning that some features will be disabled
# due to some failure.

MACRO(PRINT_WARNING message action)
  MESSAGE("WARNING: ${message}.\n   ${action}.")
ENDMACRO(PRINT_WARNING)

# Returns an unquoted string.  Note that CMake will readily turn such
# strings back into lists, due to the duality of lists and
# semicolon-separated strings.  So be careful how you use it.

MACRO(LIST2STRING alist astring)
  FOREACH(elem ${${alist}})
   SET(${astring} "${${astring}} ${elem}")
  ENDFOREACH(elem)
ENDMACRO(LIST2STRING)

# -------------------------------------------------------------
# Which modules to build?
# -------------------------------------------------------------

# For each SUNDIALS solver available (i.e. for which we have the
# sources), give the user the option of enabling/disabling it.

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")
  OPTION(BUILD_CVODE  "Build the CVODE library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")
  SET(BUILD_CVODE OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvode")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")
  OPTION(BUILD_CVODES  "Build the CVODES library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")
  SET(BUILD_CVODES OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cvodes")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")
  OPTION(BUILD_IDA  "Build the IDA library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")
  SET(BUILD_IDA OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/ida")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")
  OPTION(BUILD_IDAS  "Build the IDAS library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")
  SET(BUILD_IDAS OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/idas")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")
  OPTION(BUILD_KINSOL  "Build the KINSOL library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")
  SET(BUILD_KINSOL OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/kinsol")

IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")
  OPTION(BUILD_CPODES  "Build the CPODES library"  ON)
ELSE(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")
  SET(BUILD_CPODES OFF)
ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/src/cpodes")

# -------------------------------------------------------------
# Other configuration options
# -------------------------------------------------------------

# Option that allows users to build static and/or shared libraries
# ----------------------------------------------------------------

OPTION(BUILD_STATIC_LIBS "Build static libraries" ON)
OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Under Windows, prepare substitution variable SUNDIALS_EXPORT for sundials_config.h
IF(BUILD_SHARED_LIBRARY AND WIN32)

  # When building shared SUNDIALS libraries under Windows, use
  #      #define SUNDIALS_EXPORT __declspec(dllexport)
  # When linking to shared SUNDIALS libraries under Windows, use
  #      #define SUNDIALS_EXPORT __declspec(dllimport)
  SET(SUNDIALS_EXPORT
    "#ifdef BUILD_SUNDIALS_LIBRARY
#define SUNDIALS_EXPORT __declspec(dllexport)
#else
#define SUNDIALS_EXPORT __declspec(dllimport)
#endif")

ELSE(BUILD_SHARED_LIBRARY AND WIN32)

  # In all other cases (other platforms or static libraries
  # under Windows), the SUNDIALS_EXPORT macro is empty
  SET(SUNDIALS_EXPORT "#define SUNDIALS_EXPORT")

ENDIF(BUILD_SHARED_LIBRARY AND WIN32)

# Make sure we build at least one type of libraries
IF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  PRINT_WARNING("Both static and shared library generation were disabled"
    "Building static libraries was re-enabled")
  SET(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
ENDIF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)

# Option to specify precision
# ---------------------------

SET(SUNDIALS_PRECISION "double" CACHE STRING "double, single or extended")

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
STRING(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
SET(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

# Option to use the generic math libraries (UNIX only)
# ----------------------------------------------------

IF(UNIX)
  OPTION(USE_GENERIC_MATH "Use generic (std-c) math libraries" ON)
  IF(USE_GENERIC_MATH)
    # executables will be linked against -lm 
    SET(EXTRA_LINK_LIBS -lm)
    # prepare substitution variable GENERIC_MATH_LIB for sundials_config.h
    SET(GENERIC_MATH_LIB "#define SUNDIALS_USE_GENERIC_MATH 1")
  ENDIF(USE_GENERIC_MATH)
ENDIF(UNIX)

# -------------------------------------------------------------
# Variable initializations
# -------------------------------------------------------------

# Each of the variables below (of type FEATURE_STATUS) can have 
# one of the following values:
#  none         initial value indicating that the feature is
#               not needed for the current configuration.
#  needed       the feature is required but not all prerequisites
#               have been tested yet.
#  available    the feature is required and available (all
#               prerequisites are also available)
#  unavailable  the feature is required but it is not available
#               (at least one prerequisite not available)

# Status of FCMIX support
SET(FCMIX_STATUS "none")

# Status of Lapack support
SET(LAPACK_STATUS "none")

# Status of MPI C
SET(MPIC_STATUS "none")

# Status of MPI F77
SET(MPIF_STATUS "none")

# Status of F77 compiler
SET(F77_STATUS "none")

# Status of F77 name-mangling scheme
SET(SCHEME_STATUS "none")

# -------------------------------------------------------------
# Enable Fortran support?
# -------------------------------------------------------------

# FCMIX support is an option only if at least one solver that 
# provides such an interface is built.

IF(BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)
  SHOW_VARIABLE(FCMIX_ENABLE BOOL "Enable Fortran-C support" OFF)
  IF(FCMIX_ENABLE)
    SET(FCMIX_STATUS "needed")
  ENDIF(FCMIX_ENABLE)
ELSE(BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)
  HIDE_VARIABLE(FCMIX_ENABLE)
ENDIF(BUILD_CVODE OR BUILD_IDA OR BUILD_KINSOL)

# -------------------------------------------------------------
# Enable BLAS/LAPACK support?
# -------------------------------------------------------------

OPTION(LAPACK_ENABLE "Enable Lapack support" OFF)

IF(LAPACK_ENABLE)
  SET(LAPACK_STATUS "needed")
ENDIF(LAPACK_ENABLE)

# -------------------------------------------------------------
# Enable MPI support?
# -------------------------------------------------------------

OPTION(MPI_ENABLE "Enable MPI support" OFF)

IF(MPI_ENABLE)
  SET(MPIC_STATUS "needed")
  IF(FCMIX_ENABLE)
    SET(MPIF_STATUS "needed")
  ENDIF(FCMIX_ENABLE)
ENDIF(MPI_ENABLE)

# -------------------------------------------------------------
# Enable examples?
# -------------------------------------------------------------

OPTION(EXAMPLES_ENABLE "Build the SUNDIALS examples" OFF)

IF(EXAMPLES_ENABLE)

  # If examples are enabled, set different options

  # The examples will be linked with the library corresponding 
  # to the build type
  IF(BUILD_SHARED_LIBS)
    SET(LINK_LIBRARY_TYPE "shared")
  ELSE(BUILD_SHARED_LIBS)
    SET(LINK_LIBRARY_TYPE "static")
  ENDIF(BUILD_SHARED_LIBS)

  # However, if both shared and static SUNDIALS libraries are built,
  # the default is to link the examples against the shared libraries.
  # the user can change this thourgh the option EXAMPLES_USE_STATIC_LIBS  
  IF(BUILD_STATIC_LIBS AND BUILD_SHARED_LIBS)
    SHOW_VARIABLE(EXAMPLES_USE_STATIC_LIBS BOOL
      "Link examples with the shared SUNDIALS libraries" OFF)
    IF(EXAMPLES_USE_STATIC_LIBS)
      SET(LINK_LIBRARY_TYPE "static")
    ENDIF(EXAMPLES_USE_STATIC_LIBS)
  ELSE(BUILD_STATIC_LIBS AND BUILD_SHARED_LIBS)
    HIDE_VARIABLE(EXAMPLES_USE_STATIC_LIBS)
  ENDIF(BUILD_STATIC_LIBS AND BUILD_SHARED_LIBS)

  # Check if example files are to be exported
  SHOW_VARIABLE(EXAMPLES_INSTALL BOOL
    "Install example files" ON)

  # If examples are to be exported, check where we should install them.
  IF(EXAMPLES_INSTALL)

    SHOW_VARIABLE(EXAMPLES_INSTALL_PATH STRING
      "Output directory for installing example files" "${CMAKE_INSTALL_PREFIX}/examples")
    
    IF(NOT EXAMPLES_INSTALL_PATH)
      PRINT_WARNING("The example installation path is empty"
        "Example installation path was reset to its default value")
      SET(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/examples" CACHE STRING
        "Output directory for installing example files" FORCE)
    ENDIF(NOT EXAMPLES_INSTALL_PATH)

  ELSE(EXAMPLES_INSTALL)

    HIDE_VARIABLE(EXAMPLES_INSTALL_PATH)

  ENDIF(EXAMPLES_INSTALL)

  # On UNIX systems, if we install examples, check if we should also
  # generate and install Makefiles.
  IF(UNIX AND EXAMPLES_INSTALL)
    SHOW_VARIABLE(EXAMPLES_GENERATE_MAKEFILES BOOL
      "Create Makefiles for building the examples" ON)
  ELSE(UNIX AND EXAMPLES_INSTALL)
    HIDE_VARIABLE(EXAMPLES_GENERATE_MAKEFILES)
  ENDIF(UNIX AND EXAMPLES_INSTALL)

ELSE(EXAMPLES_ENABLE)

  # If examples are disabled, hide all options related to
  # building and installing the SUNDIALS examples

  HIDE_VARIABLE(EXAMPLES_INSTALL)
  HIDE_VARIABLE(EXAMPLES_INSTALL_PATH)
  HIDE_VARIABLE(EXAMPLES_USE_STATIC_LIBS)
  HIDE_VARIABLE(EXAMPLES_GENERATE_MAKEFILES)
  
ENDIF(EXAMPLES_ENABLE)


# -------------------------------------------------------------
# Save compiler flags corresponding to the current build type
# -------------------------------------------------------------

# Make sure that the following tests use the C and Fortran flags corresponding
# to the current build type. These flags are stored in the variables TMP_C_FLAGS 
# and TMP_Fortran_FLAGS, respectively, and are used in the generated CMakeLists files.
IF(NOT CMAKE_BUILD_TYPE)
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
ENDIF(NOT CMAKE_BUILD_TYPE)
IF(CMAKE_BUILD_TYPE MATCHES "Default")
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
ENDIF(CMAKE_BUILD_TYPE MATCHES "Default")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS_RELEASE})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
ENDIF(CMAKE_BUILD_TYPE MATCHES "Release")
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS_DEBUG})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")
IF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS_RELWITHDEBINFO})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})
ENDIF(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
IF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
  SET(TMP_C_FLAGS       ${CMAKE_C_FLAGS_MINSIZE})
  SET(TMP_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_MINSIZE})
ENDIF(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")

# -------------------------------------------------------------
# Test Fortran compiler with current flags
# -------------------------------------------------------------

# A Fortran compiler is "needed" if FCMIX_STATUS is "needed"
# or LAPCK_STATUS is "needed" (in the later case only because
# the F77 name-manling scheme is required and we'd like to
# be able to determine it automatically, in which case we
# need a working Fortran compiler).
# Therefore, if a working compiler is not found, we set
# FCMIX_STATUS to "unavailable" but we do not change 
# LAPACK_STATUS.

IF(FCMIX_STATUS STREQUAL "needed" OR LAPACK_STATUS STREQUAL "needed")
  SET(F77_STATUS "needed")
ENDIF(FCMIX_STATUS STREQUAL "needed" OR LAPACK_STATUS STREQUAL "needed")

IF(F77_STATUS STREQUAL "needed")

  # Test if we can compile and link a simple Fortran program.
  # We create first a library flib defining two functions (mysub and my_sub)
  # which will also be used later to determine the Fortran name-mangling scheme.
  # Nest we attempt to create an executable which calles these two functions.

  MESSAGE(STATUS "Determining if the Fortran compiler works...")

  # Enable Fortran support
  ENABLE_LANGUAGE(Fortran)
  MARK_AS_ADVANCED(CLEAR 
    CMAKE_Fortran_COMPILER
    CMAKE_Fortran_FLAGS
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_MINSIZEREL
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_RELWITHDEB)

  # Create the FortranTest directory
  SET(FortranTest_DIR ${PROJECT_BINARY_DIR}/FortranTest)
  FILE(MAKE_DIRECTORY ${FortranTest_DIR})

  # Create a CMakeLists.txt file which will generate the "flib" library
  # and an executable "ftest"
  FILE(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "PROJECT(FortranTest Fortran)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
    "ADD_LIBRARY(flib flib.f)\n"
    "ADD_EXECUTABLE(ftest ftest.f)\n"
    "TARGET_LINK_LIBRARIES(ftest flib)\n")
  
  # Create the Fortran source flib.f which defines two subroutines, "mysub" and "my_sub"
  FILE(WRITE ${FortranTest_DIR}/flib.f
    "        SUBROUTINE mysub\n"
    "        RETURN\n"
    "        END\n"
    "        SUBROUTINE my_sub\n"
    "        RETURN\n"
    "        END\n")

  # Create the Fortran source ftest.f which calls "mysub" and "my_sub"
  FILE(WRITE ${FortranTest_DIR}/ftest.f
    "        PROGRAM ftest\n"
    "        CALL mysub()\n"
    "        CALL my_sub()\n"
    "        END\n")
  
  # Use TRY_COMPILE to make the targets "flib" and "ftest"
  TRY_COMPILE(FTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
    ftest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)

  # Change F77_STATUS based on test result. If test fails this
  # automatically results in disabling FCMIX and MPI-F support.
  IF(FTEST_OK)
    SET(F77_STATUS "available")
  ELSE(FTEST_OK)
    SET(F77_STATUS "unavailable")
    IF(FCMIX_STATUS STREQUAL "needed")
      PRINT_WARNING("Unable to compile and link a simple Fortran program"
        "FCMIX support will not be provided")
      SET(FCMIX_STATUS "unavailable")
    ENDIF(FCMIX_STATUS STREQUAL "needed")
    IF(MPIF_STATUS STREQUAL "needed")
      PRINT_WARNING("Unable to compile and link a simple Fortran program"
        "MPI Fortran support will not be provided")
      SET(MPIF_STATUS "unavailable")
    ENDIF(MPIF_STATUS STREQUAL "needed")
  ENDIF(FTEST_OK)

ENDIF(F77_STATUS STREQUAL "needed")

# -------------------------------------------------------------
# Set the Fortran name-mangling scheme if needed 
# -------------------------------------------------------------

# The F77 name-mangling scheme is "needed" if FCMIX_STATUS is
# "needed" or if LAPACK_STATUS is "needed".

IF(FCMIX_STATUS STREQUAL "needed" OR LAPACK_STATUS STREQUAL "needed")
  SET(SCHEME_STATUS "needed")
ENDIF(FCMIX_STATUS STREQUAL "needed" OR LAPACK_STATUS STREQUAL "needed")

# NOTES:
#
# In general, names of symbols with and without underscore may be mangled 
# differently (e.g. g77 mangles mysub to mysub_ and my_sub to my_sub__),
# we have to consider both cases.
#
# Once the name-mangling schemes are determined, we use them to define two C
# preprocessor macros, F77_FUNC and F77_FUNC_, corresponding to the two cases:
# symbols with names not containing underscores and symbols with names containing 
# underscores. For example, if using g77 the definitions of these two macros will be:
#     #define F77_FUNC(name,NAME) name ## _
#     #define F77_FUNC_(name,NAME) name ## __
# The appropriate #define lines are stored in the variables DEFINE_F77_FUNC and
# DEFINE_F77_FUNC_, respectively, and can be used to generate a configuration header
# file (using CONFIGURE_FILE).
# 
# Since all FCMIX interface functions in SUNDIALS use function names NOT containing 
# underscores, in case we cannot determine the scheme automatically, we simply ignore 
# the case of symbol names containing underscores (to minimize user input).

# If the name-mangling scheme is required AND a working Fortran compiler is available, 
# try to automatically determine the F77 name-mangling scheme.

IF(SCHEME_STATUS STREQUAL "needed" AND F77_STATUS STREQUAL "available")

  # (A) Automatic detection of Fortran name-mangling scheme
  #     ---------------------------------------------------
  #
  #  1) create a library from a Fortran source file which defines a function "mysub"
  #  2) attempt to link with this library a C source file which calls the "mysub"
  #     function using various possible schemes (6 different schemes, corresponding
  #     to all combinations lower/upper case and none/one/two underscores)
  #  3) define the name-mangling scheme based on the test that was successful.

  MESSAGE(STATUS "Determining the Fortran name-mangling scheme...")

  SET(HAVE_SCHEME_NO_UNDERSCORES FALSE)
  SET(HAVE_SCHEME_WITH_UNDERSCORES FALSE)

  # Initialize the name-mangling schemes for symbol names with and without underscores
  # If successful indetermining the schemes, SCHEME_NO_UNDERSCORE will be set to the
  # mangled name of "mysub" and SCHEME_WITH_UNDERSORES will be set to the mangled name
  # of "my_sub".
  SET(SCHEME_NO_UNDERSCORES "")
  SET(SCHEME_WITH_UNDERSCORES "")

  # CASE 1: symbol names WITHOUT undersores
  # ---------------------------------------
  
  # Overwrite CMakeLists.txt with one which will generate the "ctest1" executable
  FILE(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "PROJECT(FortranTest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
    "ADD_EXECUTABLE(ctest1 ctest1.c)\n"
    "FIND_LIBRARY(FLIB flib ${FortranTest_DIR})\n"
    "TARGET_LINK_LIBRARIES(ctest1 \${FLIB})\n")
  
  # Define the list "options" of all possible schemes that we want to consider
  # Get its length and initialize the counter "iopt" to zero
  SET(options mysub mysub_ mysub__ MYSUB MYSUB_ MYSUB__)
  LIST(LENGTH options imax)
  SET(iopt 0)
  
  # We will attempt to sucessfully generate the "ctest" executable as long as
  # there still are entries in the "options" list
  WHILE(${iopt} LESS ${imax})   
    
    # Get the current list entry (current scheme)
    LIST(GET options ${iopt} opt)
    
    # Generate C source which calls the "mysub" function using the current scheme
    FILE(WRITE ${FortranTest_DIR}/ctest1.c "int main(){${opt}();return(0);}\n")
    
    # Use TRY_COMPILE to make the "ctest1" executable from the current C source
    # and linking to the previously created "flib" library.
    TRY_COMPILE(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
      ctest1 OUTPUT_VARIABLE MY_OUTPUT)
    
    # To ensure we do not use stuff from the previous attempts, 
    # we must remove the CMakeFiles directory.
    FILE(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
    
    # Test if we successfully created the "ctest" executable.
    # If yes, save the current scheme, and set the counter "iopt" to "imax" 
    # so that we exit the while loop.
    # Otherwise, increment the counter "iopt" and go back in the while loop.
    IF(CTEST_OK)
      SET(SCHEME_NO_UNDERSCORES ${opt})
      SET(iopt ${imax})
    ELSE(CTEST_OK)
      MATH(EXPR iopt ${iopt}+1)
    ENDIF(CTEST_OK)
    
  ENDWHILE(${iopt} LESS ${imax})   

  # If the name-mangling scheme of symbol names not containing underscores
  # was successfully determined (i.e. SCHEME_NO_UNDERSCORES is non-empty, 
  # set the appropriate C preprocessor macro

  IF(SCHEME_NO_UNDERSCORES)

    SET(SCHEME_STATUS "available")

    IF(SCHEME_NO_UNDERSCORES MATCHES "mysub")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "mysub")
    
    IF(SCHEME_NO_UNDERSCORES MATCHES "mysub_")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name ## _")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "mysub_")      
    
    IF(SCHEME_NO_UNDERSCORES MATCHES "mysub__")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name ## __")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "mysub__")
    
    IF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB")
    
    IF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB_")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME ## _")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB_")
    
    IF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB__")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME ## __")
    ENDIF(SCHEME_NO_UNDERSCORES MATCHES "MYSUB__")
    
    MESSAGE(STATUS 
      "name-mangling scheme for symbol names without underscores:\n"
      "   mysub  ->  ${SCHEME_NO_UNDERSCORES}\n"
      "   ${DEFINE_F77_FUNC}")

  ELSE(SCHEME_NO_UNDERSCORES)

    SET(DEFINE_F77_FUNC "")

  ENDIF(SCHEME_NO_UNDERSCORES)

  # CASE 2: symbol names WITH undersores
  # ------------------------------------
  
  FILE(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "PROJECT(FortranTest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
    "ADD_EXECUTABLE(ctest2 ctest2.c)\n"
    "FIND_LIBRARY(FLIB flib ${FortranTest_DIR})\n"
    "TARGET_LINK_LIBRARIES(ctest2 \${FLIB})\n")
  
  SET(options my_sub my_sub_ my_sub__ MY_SUB MY_SUB_ MY_SUB__)
  LIST(LENGTH options imax)
  SET(iopt 0)
  
  WHILE(${iopt} LESS ${imax})   
    LIST(GET options ${iopt} opt)
    FILE(WRITE ${FortranTest_DIR}/ctest2.c "int main(){${opt}();return(0);}\n")
    TRY_COMPILE(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
      ctest2 OUTPUT_VARIABLE MY_OUTPUT)
    FILE(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
    IF(CTEST_OK)
      SET(SCHEME_WITH_UNDERSCORES ${opt})
      SET(iopt ${imax})
    ELSE(CTEST_OK)
      MATH(EXPR iopt ${iopt}+1)
    ENDIF(CTEST_OK)
  ENDWHILE(${iopt} LESS ${imax})   
    
  # If the name-mangling scheme of symbol names containing underscores
  # was successfully determined (i.e. SCHEME_WITH_UNDERSCORES is non-empty, 
  # set the appropriate C preprocessor macro

  IF(SCHEME_WITH_UNDERSCORES)
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) name")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub")
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub_")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) name ## _")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub_")      
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub__")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) name ## __")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "my_sub__")
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) NAME")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB")
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB_")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) NAME ## _")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB_")
    
    IF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB__")
      SET(DEFINE_F77_FUNC_ "#define F77_FUNC_(name,NAME) NAME ## __")
    ENDIF(SCHEME_WITH_UNDERSCORES MATCHES "MY_SUB__")
    
    MESSAGE(STATUS 
      "name-mangling scheme for symbol names with underscores:\n"
      "   my_sub ->  ${SCHEME_WITH_UNDERSCORES}\n"
      "   ${DEFINE_F77_FUNC_}")

  ELSE(SCHEME_WITH_UNDERSCORES)

    SET(DEFINE_F77_FUNC_ "")

  ENDIF(SCHEME_WITH_UNDERSCORES)

ENDIF(SCHEME_STATUS STREQUAL "needed" AND F77_STATUS STREQUAL "available")


# If the name-mangling scheme is still "needed" at this stage, it means that
# either a working Fortran compiler was not found or that automatic detection
# of the scheme failed. In this case ask the user to provide it by specifying 
# the case and number of appended underscores.

IF(SCHEME_STATUS STREQUAL "needed")

  # (B) User-provided Fortran name-mangling scheme
  #     ------------------------------------------
  # 
  # Ask the user to specify the case and number of appended underscores
  # corresponding to the Fortran name-mangling scheme of symbol names 
  # that do not themselves contain underscores (recall that this is all
  # we really need for the SUNDIALS FCMIX interfaces).
  # Note: the default scheme is lower case - one underscore

  # Specify the case for the Fortran name-mangling scheme
  SHOW_VARIABLE(F77_FUNC_CASE STRING
    "case of Fortran function names (lower/upper)"
    "lower")
  # Specify the number of appended underscores for the Fortran name-mangling scheme
  SHOW_VARIABLE(F77_FUNC_UNDERSCORES STRING 
    "number of underscores appended to Fortran function names"
    "one")
  
  # Based on the given case and number of underscores,
  # set the C preprocessor macro definition

  IF(${F77_FUNC_CASE} MATCHES "lower")
    
    IF(${F77_FUNC_UNDERSCORES} MATCHES "none")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "none")
    IF(${F77_FUNC_UNDERSCORES} MATCHES "one")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name ## _")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "one")
    IF(${F77_FUNC_UNDERSCORES} MATCHES "two")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) name ## __")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "two")
    
  ELSE(${F77_FUNC_CASE} MATCHES "lower")

    IF(${F77_FUNC_UNDERSCORES} MATCHES "none")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "none")
    IF(${F77_FUNC_UNDERSCORES} MATCHES "one")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME ## _")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "one")
    IF(${F77_FUNC_UNDERSCORES} MATCHES "two")
      SET(DEFINE_F77_FUNC "#define F77_FUNC(name,NAME) NAME ## __")
    ENDIF(${F77_FUNC_UNDERSCORES} MATCHES "two")
    
  ENDIF(${F77_FUNC_CASE} MATCHES "lower")

  # Test the user-provided scheme 
  # (use the F77_FUNC macro defined in DEFINE_F77_FUNC)

  FILE(WRITE ${FortranTest_DIR}/CMakeLists.txt
    "PROJECT(FortranTest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
    "ADD_EXECUTABLE(ctest3 ctest3.c)\n"
    "FIND_LIBRARY(FLIB flib ${FortranTest_DIR})\n"
    "TARGET_LINK_LIBRARIES(ctest3 \${FLIB})\n")
  
  FILE(WRITE ${FortranTest_DIR}/ctest3.c
    "${DEFINE_F77_FUNC}\n"
    "#define sub F77_FUNC(mysub, MYSUB)\n"
    "int main(){\n"
    "sub();\n"
    "return(0);\n"
    "}\n")
  
  TRY_COMPILE(CTEST_OK ${FortranTest_DIR} ${FortranTest_DIR}
    ctest3 OUTPUT_VARIABLE MY_OUTPUT)
  
  FILE(REMOVE_RECURSE ${FortranTest_DIR}/CMakeFiles)
  
  IF(CTEST_OK)
    SET(SCHEME_STATUS "available")
  ELSE(CTEST_OK)
    SET(SCHEME_STATUS "unavailable")
    SET(DEFINE_F77_FUNC "")
  ENDIF(CTEST_OK)

  # Since the SUNDIALS FCMIX modules do not use symbol names containing
  # underscores, simply set DEFINE_F77_FUNC_ to the empty string
  SET(DEFINE_F77_FUNC_ "")

ELSE(SCHEME_STATUS STREQUAL "needed")

  # The Fortran name-mangling scheme is either not needed
  # or else it was automatically determined.

  # Hide options which are not needed.
  HIDE_VARIABLE(F77_FUNC_CASE)
  HIDE_VARIABLE(F77_FUNC_UNDERSCORES)

ENDIF(SCHEME_STATUS STREQUAL "needed")

# If FCMIX_STATUS is "needed" change its value to "available 
# if the F77 name-mangling scheme is available

IF(FCMIX_STATUS STREQUAL "needed" AND SCHEME_STATUS STREQUAL "available")
  SET(FCMIX_STATUS "available")
ENDIF(FCMIX_STATUS STREQUAL "needed" AND SCHEME_STATUS STREQUAL "available")

# If the F77 name-mangling scheme is "unavailable" at this point, this
# means that although required, we were unable to determine it automatically
# and the user provided scheme is incorrect.
# This automatically results in disabling FCMIX, MPI-F, and LAPACK support.

IF(SCHEME_STATUS STREQUAL "unavailable")

  IF(FCMIX_STATUS STREQUAL "needed")
    PRINT_WARNING("Unable to determine the Fortran name-mangling scheme"
      "FCMIX support will not be provided")
    SET(FCMIX_STATUS "unavailable")
  ENDIF(FCMIX_STATUS STREQUAL "needed")

  IF(MPIF_STATUS STREQUAL "needed")
    PRINT_WARNING("Unable to determine the Fortran name-mangling scheme"
      "MPI Fortran support will not be provided")
    SET(MPIF_STATUS "unavailable")
  ENDIF(MPIF_STATUS STREQUAL "needed")

  IF(LAPACK_STATUS STREQUAL "needed")
    PRINT_WARNING("Unable to determine the Fortran name-mangling scheme"
      "LAPACK support will not be provided")
    SET(LAPACK_STATUS "unavailable")
  ENDIF(LAPACK_STATUS STREQUAL "needed")

ENDIF(SCHEME_STATUS STREQUAL "unavailable")

# Copy DEFINE_F77_FUNC and DEFINE_F77_FUNC_ in the substitution variables 
# F77_MANGLE_MACRO1 and F77_MANGLE_MACRO2, as required by the template 
# sundials_config.in (if the F77 name-mangling scheme is not needed or 
# unavailable, these two variables are empty strings)

SET(F77_MANGLE_MACRO1 ${DEFINE_F77_FUNC})
SET(F77_MANGLE_MACRO2 ${DEFINE_F77_FUNC_})

# -------------------------------------------------------------
# Find and test the Lapack libraries
# -------------------------------------------------------------

# If LAPACK_STATUS is "needed", first try to find the appropriate
# libraries and linker flags needed to link against them.

IF(LAPACK_STATUS STREQUAL "needed")

  MESSAGE(STATUS "Checking for Lapack libraries...")

  IF(F77_STATUS STREQUAL "available")
    # Look for libraries only if none are defined
    IF(NOT LAPACK_LIBRARIES)
      INCLUDE(FindLAPACK)
    ENDIF(NOT LAPACK_LIBRARIES)
  ENDIF(F77_STATUS STREQUAL "available")

  # If using a GNU C compiler, it is quite likely we'll want LAPACK_LINKER_FLAGS
  # to include -lg2c (if not already present)
  IF(CMAKE_COMPILER_IS_GNUCC)
    IF(NOT LAPACK_LINKER_FLAGS MATCHES "g2c")
      SET(LAPACK_LINKER_FLAGS "${LAPACK_LINKER_FLAGS} -lg2c")
    ENDIF(NOT LAPACK_LINKER_FLAGS MATCHES "g2c")
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  # Expose to the user (and cache) LAPACK_LIBRARIES and LAPACK_LINKER_FLAGS
  SHOW_VARIABLE(LAPACK_LIBRARIES STRING
    "Lapack libraries" "${LAPACK_LIBRARIES}")
  SHOW_VARIABLE(LAPACK_LINKER_FLAGS STRING
    "Lapack required linker flags" "${LAPACK_LINKER_FLAGS}")

  # If the LAPACK_LIBRARIES were not found, set LAPACK_STATUS to "unavailable".
  IF(NOT LAPACK_LIBRARIES)
    PRINT_WARNING("Cannot find LAPACK libraries"
      "LAPACK support will not be provided")
    SET(LAPACK_STATUS "unavailable")
  ENDIF(NOT LAPACK_LIBRARIES)

ENDIF(LAPACK_STATUS STREQUAL "needed")

# If LAPACK_STATUS is still "needed" at this point then test libraries

IF(LAPACK_STATUS STREQUAL "needed")

  MESSAGE(STATUS "Checking if Lapack works...")

  # Create the LapackTest directory
  SET(LapackTest_DIR ${PROJECT_BINARY_DIR}/LapackTest)
  FILE(MAKE_DIRECTORY ${LapackTest_DIR})

  # Create a CMakeLists.txt file 
  FILE(WRITE ${LapackTest_DIR}/CMakeLists.txt
    "PROJECT(LapackTest C)\n"
    "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
    "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
    "SET(CMAKE_EXE_LINKER_FLAGS \"\${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}\")\n"
    "ADD_EXECUTABLE(ltest ltest.c)\n"
    "TARGET_LINK_LIBRARIES(ltest ${LAPACK_LIBRARIES})\n")

  # Create a C source file which calls a Blas function (dcopy) and an Lapack function (dgetrf)
  FILE(WRITE ${LapackTest_DIR}/ltest.c
    "${DEFINE_F77_FUNC}\n"
    "#define dcopy_f77 F77_FUNC(dcopy, DCOPY)\n"
    "#define dgetrf_f77 F77_FUNC(dgetrf, DGETRF)\n"
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
  TRY_COMPILE(LTEST_OK ${LapackTest_DIR} ${LapackTest_DIR}
    ltest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${LapackTest_DIR}/CMakeFiles)

  # Set LAPACK_STATUS based on test results
  IF(LTEST_OK)
    SET(LAPACK_STATUS "available")
  ELSE(LTEST_OK)
    PRINT_WARNING("Unable to link a C code to the Lapack libraries"
      "LAPACK support will not be provided")
    SET(LAPACK_STATUS "unavailable")
  ENDIF(LTEST_OK)

ENDIF(LAPACK_STATUS STREQUAL "needed")

# -------------------------------------------------------------
# Decide how to compile MPI codes.
# -------------------------------------------------------------

# If parallel support is "needed", check if we will be using
# MPI compiler scripts. If so, find the mpicc executable (and
# make the MPI_USE_MPISCRIPTS variable visible) and hide all variables
# relevant to the case in which we use the serial compiler with
# MPI libraries to build parallel libraries/programs. Also, if FCMIX
# is "needed", look for mpif77.
# If we do not use MPI compiler scripts, hide the MPI_MPICC variable
# and attempt to determine the location of MPI header files and the
# MPI libraries.
#
# If parallel support is disabled, hide all MPI-related variables.

IF(MPIC_STATUS STREQUAL "needed")

  MESSAGE(STATUS "Checking how to build parallel programs...")

  SHOW_VARIABLE(MPI_USE_MPISCRIPTS BOOL "Use MPI script" ON)
  
  IF(MPI_USE_MPISCRIPTS)

    # We will be using MPI compiler scripts.
    
    HIDE_VARIABLE(MPI_INCLUDE_PATH)
    HIDE_VARIABLE(MPI_LIBRARIES)
    
    FIND_PROGRAM(MPI_MPICC NAMES mpicc DOC "mpicc program")
    SHOW_VARIABLE(MPI_MPICC STRING "mpicc program" "${MPI_MPICC}")

    IF(NOT MPI_MPICC)
      PRINT_WARNING("Cannot find the MPI-C compiler script"
        "MPI support will not be provided")
      SET(MPIC_STATUS "unavailable")
      IF(MPIF_STATUS STREQUAL "needed")
        SET(MPIF_STATUS "unavailable")
      ENDIF(MPIF_STATUS STREQUAL "needed")
    ENDIF(NOT MPI_MPICC)

    IF(MPIF_STATUS STREQUAL "needed")
      FIND_PROGRAM(MPI_MPIF77 NAMES mpif77 DOC "mpif77 program")
      SHOW_VARIABLE(MPI_MPIF77 STRING "mpif77 program" "${MPI_MPIF77}")
      IF(NOT MPI_MPIF77)
        PRINT_WARNING("Cannot find MPI-F77 compiler script"
          "MPI Fortran support will not be provided")
        SET(MPIF_STATUS "unavailable")
      ENDIF(NOT MPI_MPIF77)
    ELSE(MPIF_STATUS STREQUAL "needed")
      HIDE_VARIABLE(MPI_MPIF77)
    ENDIF(MPIF_STATUS STREQUAL "needed")

  ELSE(MPI_USE_MPISCRIPTS)
    
    HIDE_VARIABLE(MPI_MPICC)
    HIDE_VARIABLE(MPI_MPIF77)
    
    INCLUDE(FindMPI)
    SHOW_VARIABLE(MPI_INCLUDE_PATH STRING "path to MPI header files" "${MPI_INCLUDE_PATH}")
    SHOW_VARIABLE(MPI_LIBRARIES STRING "MPI libraries" "${MPI_LIBRARIES}")
    
    IF(NOT MPI_LIBRARIES)
      PRINT_WARNING("Cannot find MPI libraries"
        "MPI support will not be provided")
      SET(MPIC_STATUS "unavailable")
      IF(MPIF_STATUS STREQUAL "needed")
        SET(MPIF_STATUS "unavailable")
      ENDIF(MPIF_STATUS STREQUAL "needed")
    ENDIF(NOT MPI_LIBRARIES)

  ENDIF(MPI_USE_MPISCRIPTS)

ELSE(MPIC_STATUS STREQUAL "needed")
  
  HIDE_VARIABLE(MPI_USE_MPISCRIPTS)
  HIDE_VARIABLE(MPI_MPICC)
  HIDE_VARIABLE(MPI_MPIF77)
  HIDE_VARIABLE(MPI_INCLUDE_PATH)
  HIDE_VARIABLE(MPI_LIBRARIES)
  
ENDIF(MPIC_STATUS STREQUAL "needed")

# -------------------------------------------------------------
# Test MPI
# -------------------------------------------------------------

# If MPIC_STATUS is still "needed" at this point, test C parallel programs.

IF(MPIC_STATUS STREQUAL "needed")

  MESSAGE(STATUS "Checking if MPI-C works...")

  # Create the MPITest directory
  SET(MPITest_DIR ${PROJECT_BINARY_DIR}/MPITest)
  FILE(MAKE_DIRECTORY ${MPITest_DIR})

  # Create a CMakeLists.txt file which will generate the "mpictest" executable
  IF(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPITest C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
      "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
      "ADD_EXECUTABLE(mpictest mpictest.c)\n")
  ELSE(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPITest C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpictest mpictest.c)\n"
      "TARGET_LINK_LIBRARIES(mpictest ${MPI_LIBRARIES})\n")
  ENDIF(MPI_USE_MPISCRIPTS)

  # Create a simple C source which only calls the MPI_Init and MPI_Finalize functions
  FILE(WRITE ${MPITest_DIR}/mpictest.c
    "#include <mpi.h>\n"
    "int main(){\n"
    "int c;\n"
    "char **v;\n"
    "MPI_Init(&c, &v);\n"
    "MPI_Finalize();\n"
    "return(0);\n"
    "}\n")

  # Use TRY_COMPILE to make the target "mpictest"
  TRY_COMPILE(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpitest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

  # If the test succeeds, set MPIC_STATUS to "available"
  # If the test failed, disable all MPI support (both C and Fortran)
  IF(MPITEST_OK)
    SET(MPIC_STATUS "available")
  ELSE(MPITEST_OK)
    PRINT_WARNING("Unable to compile a simple MPI C program"
      "MPI support will not be provided")
    SET(MPIC_STATUS "unavailable")
    IF(MPIF_STATUS STREQUAL "needed")
      SET(MPIF_STATUS "unavailable")
    ENDIF(MPIF_STATUS STREQUAL "needed")
  ENDIF(MPITEST_OK)

ENDIF(MPIC_STATUS STREQUAL "needed")

# If MPIF_STATUS is still "needed" at this point, test Fortran parallel programs.

IF(MPIF_STATUS STREQUAL "needed")

  MESSAGE(STATUS "Checking if MPI-F77 works...")

  # Create a CMakeLists.txt file which will generate the "mpiftest" executable
  IF(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPITest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_COMPILER ${MPI_MPIF77})\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "ADD_EXECUTABLE(mpiftest mpiftest.f)\n")
  ELSE(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPITest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpiftest mpiftest.f)\n"
      "TARGET_LINK_LIBRARIES(mpiftest ${MPI_LIBRARIES})\n")
  ENDIF(MPI_USE_MPISCRIPTS)

  # Create a simple C source which only calls the MPI_Init and MPI_Finalize functions
  FILE(WRITE ${MPITest_DIR}/mpiftest.f
    "       INCLUDE \"mpif.h\"\n"
    "       INTEGER IER\n" 
    "       CALL MPI_INIT(IER)\n"
    "       CALL MPI_FINALIZE(IER)\n"
    "       STOP\n"
    "       END\n")

  # Use TRY_COMPILE to make the target "mpiftest"
  TRY_COMPILE(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpitest OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

  # If the test succeeds, set MPIF_STATUS to "available".
  # If the test failed, disable Fortran MPI support.
  IF(MPITEST_OK)
    SET(MPIF_STATUS "available")
  ELSE(MPITEST_OK)
    PRINT_WARNING("Unable to compile a simple MPI F77 program"
      "MPI Fortran support will not be provided")
    SET(MPIF_STATUS "unavailable")
  ENDIF(MPITEST_OK)

ENDIF(MPIF_STATUS STREQUAL "needed")

# -------------------------------------------------------------
# Determine if MPI provides MPI-2 support
# -------------------------------------------------------------

# This test is not a prerequisite for MPIF. Therefore, we perform
# it only if MPI_STATUS is "available".

IF(MPIF_STATUS STREQUAL "available")

  MESSAGE(STATUS "Determining whether MPI-2 support is provided...")

  # Test for MPI_Comm_f2c() function:
  #   (1) NO  : FNVECTOR_PARALLEL module will NOT allow user to specify
  #             an MPI communicator and MPI_COMM_WORLD will be used
  #   (2) YES : FNVECTOR_PARALLEL module will allow user to specify
  #             an MPI communicator

  # Create a CMakeLists.txt file which will generate the "mpi2test" executable
  IF(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPI2Test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_COMPILER ${MPI_MPICC})\n"
      "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n")
  ELSE(MPI_USE_MPISCRIPTS)
    FILE(WRITE ${MPITest_DIR}/CMakeLists.txt
      "PROJECT(MPI2Test C)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_C_FLAGS \"${TMP_C_FLAGS}\")\n"
      "INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})\n"
      "ADD_EXECUTABLE(mpi2test mpi2test.c)\n"
      "TARGET_LINK_LIBRARIES(mpi2test ${MPI_LIBRARIES})\n")
  ENDIF(MPI_USE_MPISCRIPTS)

  # Create a simple C source which calls the MPI_Comm_f2c function
  FILE(WRITE ${MPITest_DIR}/mpi2test.c
    "#include <mpi.h>\n"
    "int main(){\n"
    "int c;\n"
    "char **v;\n"
    "MPI_Comm C_comm;\n"
    "MPI_Init(&c, &v);\n"
    "C_comm = MPI_Comm_f2c((MPI_Fint) 1);\n"
    "MPI_Finalize();\n"
    "return(0);\n"
    "}\n")

  # Use TRY_COMPILE to make the target "mpi2test"
  TRY_COMPILE(MPITEST_OK ${MPITest_DIR} ${MPITest_DIR}
    mpi2test OUTPUT_VARIABLE MY_OUTPUT)

  # To ensure we do not use stuff from the previous attempts, 
  # we must remove the CMakeFiles directory.
  FILE(REMOVE_RECURSE ${MPITest_DIR}/CMakeFiles)

  # Depending on the test result, set substitution variable for sundials_config.h
  IF(MPITEST_OK)
    MESSAGE(STATUS "MPI-2 support available")
    SET(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 1")
  ELSE(MPITEST_OK)
    MESSAGE(STATUS "MPI-2 support not available")
    SET(F77_MPI_COMM_F2C "#define SUNDIALS_MPI_COMM_F2C 0")
  ENDIF(MPITEST_OK)

ELSE(MPIF_STATUS STREQUAL "available")

  # Reset substitution variable for sundials_config.h
  SET(F77_MPI_COMM_F2C "")

ENDIF(MPIF_STATUS STREQUAL "available")

# -------------------------------------------------------------
# Configure the header file sundials_config.h
# -------------------------------------------------------------

# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )

# Add the include directory in the source tree and the one in
# the binary tree (for the header file sundials_config.h)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  )

# -------------------------------------------------------------
# Based on test results, set variables controlling the
# SUNDIALS configuration.
# -------------------------------------------------------------

# FCMIX support
IF(FCMIX_STATUS STREQUAL "available")
  SET(FCMIX_OK TRUE)
ELSE(FCMIX_STATUS STREQUAL "available")
  SET(FCMIX_OK FALSE)
ENDIF(FCMIX_STATUS STREQUAL "available")

# Lapack support
IF(LAPACK_STATUS STREQUAL "available")
  SET(LAPACK_OK TRUE)
ELSE(LAPACK_STATUS STREQUAL "available")
  SET(LAPACK_OK FALSE)
ENDIF(LAPACK_STATUS STREQUAL "available")

# MPI-C support
IF(MPIC_STATUS STREQUAL "available")
  SET(MPIC_OK TRUE)
ELSE(MPIC_STATUS STREQUAL "available")
  SET(MPIC_OK FALSE)
ENDIF(MPIC_STATUS STREQUAL "available")

# MPI-F support
IF(MPIF_STATUS STREQUAL "available")
  SET(MPIF_OK TRUE)
ELSE(MPIF_STATUS STREQUAL "available")
  SET(MPIF_OK FALSE)
ENDIF(MPIF_STATUS STREQUAL "available")

# -------------------------------------------------------------

#MESSAGE("F77 status: ${F77_STATUS}")
#MESSAGE("Scheme status: ${SCHEME_STATUS}")
#MESSAGE("FCMIX status: ${FCMIX_STATUS}")
#MESSAGE("Lapack status: ${LAPACK_STATUS}")
#MESSAGE("MPI-C status: ${MPIC_STATUS}")
#MESSAGE("MPI-F status: ${MPIF_STATUS}")
#MESSAGE("mpicc: ${MPI_MPICC}")
#MESSAGE("mpif77: ${MPI_MPIF77}")

# -------------------------------------------------------------

# -------------------------------------------------------------
# Add selected modules to the build system
# -------------------------------------------------------------

# Shared components

ADD_SUBDIRECTORY(src/sundials)
ADD_SUBDIRECTORY(src/nvec_ser)
IF(MPIC_OK)
  ADD_SUBDIRECTORY(src/nvec_par)
ENDIF(MPIC_OK)

# CVODE library

IF(BUILD_CVODE)	
  ADD_SUBDIRECTORY(src/cvode)
  IF(FCMIX_OK)
    ADD_SUBDIRECTORY(src/cvode/fcmix)
  ENDIF(FCMIX_OK)
ENDIF(BUILD_CVODE)

# CVODES library

IF(BUILD_CVODES)	
  ADD_SUBDIRECTORY(src/cvodes)
ENDIF(BUILD_CVODES)

# IDA library

IF(BUILD_IDA)	
  ADD_SUBDIRECTORY(src/ida)
  IF(FCMIX_OK)
    ADD_SUBDIRECTORY(src/ida/fcmix)
  ENDIF(FCMIX_OK)
ENDIF(BUILD_IDA)

# IDAS library

IF(BUILD_IDAS)	
  ADD_SUBDIRECTORY(src/idas)
ENDIF(BUILD_IDAS)

# KINSOL library

IF(BUILD_KINSOL)	
  ADD_SUBDIRECTORY(src/kinsol)
  IF(FCMIX_OK)
    ADD_SUBDIRECTORY(src/kinsol/fcmix)
  ENDIF(FCMIX_OK)
ENDIF(BUILD_KINSOL)

# CPODES library

IF(BUILD_CPODES)	
  ADD_SUBDIRECTORY(src/cpodes)
ENDIF(BUILD_CPODES)

# -------------------------------------------------------------
# Include the subdirectories corresponding to various examples
# -------------------------------------------------------------

# If building and installing the examples is enabled, include
# the subdirectories for those examples that will be built.
# Also, if we will generate exported example Makefiles, set 
# variables needed in generating them from templates.

IF(EXAMPLES_ENABLE)

  IF(EXAMPLES_GENERATE_MAKEFILES)
    SET(SHELL "sh")
    SET(prefix "${CMAKE_INSTALL_PREFIX}")
    SET(exec_prefix "${CMAKE_INSTALL_PREFIX}")
    SET(includedir "${prefix}/include")
    SET(libdir "${exec_prefix}/lib")
    SET(CPP "${CMAKE_C_COMPILER}")
    SET(CC "${CMAKE_C_COMPILER}")
    SET(CPPFLAGS "${CMAKE_C_FLAGS_RELEASE}")
    SET(CFLAGS "${CMAKE_C_FLAGS_RELEASE}")
    SET(LDFLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
    LIST2STRING(EXTRA_LINK_LIBS LIBS)
#    SET(LIBS "${EXTRA_LINK_LIBS}")
    IF(FCMIX_OK)
      SET(F77 "${CMAKE_Fortran_COMPILER}")
      SET(F77_LNKR "${CMAKE_Fortran_COMPILER}")
      SET(FFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      SET(F77_LDFLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
      LIST2STRING(EXTRA_LINK_LIBS F77_LIBS)
#      SET(F77_LIBS "${EXTRA_LINK_LIBS}")
    ENDIF(FCMIX_OK)
    IF(LAPACK_OK)
      LIST2STRING(LAPACK_LIBRARIES BLAS_LAPACK_LIBS)
#      SET(BLAS_LAPACK_LIBS "${LAPACK_LIBRARIES}")
      SET(BLAS_LAPACK_LIBS "${LAPACK_LINKER_FLAGS} ${BLAS_LAPACK_LIBS}")
    ENDIF(LAPACK_OK)
    IF(MPIC_OK)
      IF(MPI_USE_MPISCRIPTS)
        SET(MPICC "${MPI_MPICC}")
        IF(MPIF_OK)
          SET(MPIF77 "${MPI_MPIF77}")
          SET(MPIF77_LNKR "${MPI_MPIF77}")
        ENDIF(MPIF_OK)
        SET(MPI_INC_DIR ".")
        SET(MPI_LIB_DIR ".")
        SET(MPI_LIBS "")
        SET(MPI_FLAGS "")
      ELSE(MPI_USE_MPISCRIPTS)
        SET(MPICC "${CMAKE_C_COMPILER}")
        IF(MPIF_OK)
          SET(MPIF77 "${CMAKE_Fortran_COMPILER}")
          SET(MPIF77_LNKR "${CMAKE_Fortran_COMPILER}")
        ENDIF(MPIF_OK)
        SET(MPI_INC_DIR "${MPI_INCLUDE_PATH}")
        SET(MPI_LIB_DIR ".")
        LIST2STRING(MPI_LIBRARIES MPI_LIBS)
#        SET(MPI_LIBS "${MPI_LIBRARIES}")
      ENDIF(MPI_USE_MPISCRIPTS)
    ENDIF(MPIC_OK)
  ENDIF(EXAMPLES_GENERATE_MAKEFILES)

  IF(BUILD_CVODE)
    ADD_SUBDIRECTORY(examples/cvode/serial)
    IF(FCMIX_OK)
      ADD_SUBDIRECTORY(examples/cvode/fcmix_serial)
    ENDIF(FCMIX_OK)
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/cvode/parallel)
    ENDIF(MPIC_OK)
    IF(MPIF_OK)
      ADD_SUBDIRECTORY(examples/cvode/fcmix_parallel)
    ENDIF(MPIF_OK)
  ENDIF(BUILD_CVODE)
  
  IF(BUILD_CVODES)	
    ADD_SUBDIRECTORY(examples/cvodes/serial)
    IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/test_examples/cvodes/serial")
      ADD_SUBDIRECTORY(test_examples/cvodes/serial)
    ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/test_examples/cvodes/serial")
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/cvodes/parallel)
      IF(IS_DIRECTORY "${sundials_SOURCE_DIR}/test_examples/cvodes/parallel")
        ADD_SUBDIRECTORY(test_examples/cvodes/parallel)
      ENDIF(IS_DIRECTORY "${sundials_SOURCE_DIR}/test_examples/cvodes/parallel")
    ENDIF(MPIC_OK)
  ENDIF(BUILD_CVODES)
  
  IF(BUILD_IDA)
    ADD_SUBDIRECTORY(examples/ida/serial)
    IF(FCMIX_OK)
      ADD_SUBDIRECTORY(examples/ida/fcmix_serial)
    ENDIF(FCMIX_OK)
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/ida/parallel)
    ENDIF(MPIC_OK)
    IF(MPIF_OK)
      ADD_SUBDIRECTORY(examples/ida/fcmix_parallel)
    ENDIF(MPIF_OK)
  ENDIF(BUILD_IDA)
  
  IF(BUILD_IDAS)	
    ADD_SUBDIRECTORY(examples/idas/serial)
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/idas/parallel)
    ENDIF(MPIC_OK)
  ENDIF(BUILD_IDAS)

  IF(BUILD_KINSOL)	
    ADD_SUBDIRECTORY(examples/kinsol/serial)
    IF(FCMIX_OK)
      ADD_SUBDIRECTORY(examples/kinsol/fcmix_serial)
    ENDIF(FCMIX_OK)
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/kinsol/parallel)
    ENDIF(MPIC_OK)
    IF(MPIF_OK)
      ADD_SUBDIRECTORY(examples/kinsol/fcmix_parallel)
    ENDIF(MPIF_OK)
  ENDIF(BUILD_KINSOL)

  IF(BUILD_CPODES)	
    ADD_SUBDIRECTORY(examples/cpodes/serial)
    IF(MPIC_OK)
      ADD_SUBDIRECTORY(examples/cpodes/parallel)
    ENDIF(MPIC_OK)
  ENDIF(BUILD_CPODES)
  
ENDIF(EXAMPLES_ENABLE)

#----------------------------------
# Install configuration header file
#----------------------------------

# install configured header file
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  DESTINATION include/sundials
  )

